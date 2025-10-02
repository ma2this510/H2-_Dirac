from sacred import Experiment
from sacred.observers import FileStorageObserver, MongoObserver
from dotenv import load_dotenv
from collections import deque
import subprocess
import os
import glob
import requests
import psutil
import time

load_dotenv()

ex = Experiment("H2+_Dirac_NKB")
ex.observers.append(FileStorageObserver("experiments"))
ex.observers.append(MongoObserver(url=os.environ['URI'], db_name=os.environ['DB']))

def get_git_commit(path="fortran"):
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=path
        ).decode().strip()
    except Exception:
        return None


def extract_last_eigenvalue(filename):
    try:
        with open(filename, "rb") as f:
            one_before_last = deque(f, 2)[0]

            text = one_before_last.decode("utf-8").strip()
            parts = text.split()
            return parts[-1] if parts else None
    except Exception as e:
        print(f"Failed to read last eigenvalue: {e}")
    return None


def push(message: str):
    try:
        # Get hostname using `hostname` command
        hostname = subprocess.check_output(["hostname"], text=True).strip()

        # Send notification to ntfy.sh
        response = requests.post(
            os.environ['NTFY'],
            data=message.encode("utf-8"),
            headers={
                "Title": f"From {hostname}",
                "Priority": "default"
            },
            timeout=5
        )
        response.raise_for_status()
    except (subprocess.CalledProcessError, requests.RequestException) as e:
        print(f"Notification failed: {e}")

@ex.config
def config():
    git_commit = get_git_commit()
    d = 8
    n = 20
    n_remove = 0
    Z1 = 1.0
    Z2 = 1.0
    m = 1.0
    c = 137.035999679
    R = 1.0
    ximax = 30.0
    ximin = 1.0
    epsilon = 0.0
    eta_slp = 4.0e-2
    save_step = ".false."


@ex.named_config
def dithorium():
    git_commit = get_git_commit()
    d = 8
    n = 20
    n_remove = 0
    Z1 = 90.0
    Z2 = 90.0
    m = 1.0
    c = 137.035999679
    R = 1.11e-2
    ximax = 30.0
    ximin = 1.0
    epsilon = 0.0
    eta_slp = 4.0e-2
    save_step = ".false."


@ex.automain
def run(d, n, n_remove, Z1, Z2, m, c, R, ximax, ximin, epsilon, eta_slp, save_step):
    # Write formatted input file
    with open("input.txt", "w") as f:
        f.write(f"""\
    {d},                   - d : Order of the B-Spline (order Mathematica + 1)
    {n},                  - n : Number of used B-Splines
    {n_remove},                   - n_remove : Number of B-Splines to remove at each end
    {Z1:.2f},                - Z1 : Z1 value
    {Z2:.2f},                - Z2 : Z2 value
    {m:.2f},                - m : electron mass
    {c:.9f},       - c : speed of light
    {R:.1f},                 - R : SEMI-interatomic distance
    {ximax:.1f},                - ximax : Maximum xi value
    {ximin:.1f},                 - ximin : Minimum xi value
    {epsilon:.1f},                 - epsilon : to avoid singularities
    {eta_slp:.1e},              - eta_slp : slope value to generate knots on eta
    {save_step},              - save_step : Save every matrices
""")

    # Compile and run Fortran program
    print("Compiling and running Fortran code...")
    subprocess.run(["make", "-C", "fortran", "main.out"], check=True, text=True)
    # result = subprocess.run(["./fortran/main.out"], input=open("input.txt").read(), text=True, capture_output=True)
    result = subprocess.Popen(["./fortran/main.out"], stdin=open("input.txt"),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Monitor memory and CPU usage
    pid = result.pid
    p = psutil.Process(pid)
    print(f"Monitoring process PID: {pid}")
    p.cpu_percent(interval=None) # Initialize CPU percent calculation
    try:
        while True:
            ex.log_scalar("memory_mb", p.memory_info().rss / (1024 * 1024))
            ex.log_scalar("cpu_percent", p.cpu_percent(interval=None))
            if result.poll() is not None:
                break
            time.sleep(5)
    except psutil.NoSuchProcess:
        pass

    # Save raw logs and outputs
    with open("output.log", "w") as f:
        f.write(result.stdout.read())

    ex.add_artifact("input.txt")
    ex.add_artifact("output.log")
    ex.add_artifact("eigenvalues.txt")
    ex.add_artifact("log_file")

    for file in glob.glob("*.csv"):
        ex.add_artifact(file)

    branch = subprocess.check_output(
        ["git", "-C", "fortran", "rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
    cwd = subprocess.check_output(["pwd"]).decode().strip()
    msg = f"✅ Computation finished on branch \"{branch}\" in directory \"{cwd}\""
    push(msg)

    last_eigenvalue = extract_last_eigenvalue("eigenvalues.txt")
    return {"last_eigenvalue": last_eigenvalue}
