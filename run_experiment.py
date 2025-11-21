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
ex.observers.append(MongoObserver(
    url=os.environ['URI'], db_name=os.environ['DB']))


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
    xi_slp = 7.0
    save_step = ".false."
    tot_diag = ".false."
    maxit = 10
    eig = 1.87777625e4
    compute_wf = ".false."


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
    R = 0.011111
    ximax = 20.0
    ximin = 1.0
    epsilon = 0.0
    eta_slp = 0.029
    xi_slp = 4.0e1
    save_step = ".false."
    tot_diag = ".false."
    maxit = 10
    eig = 9.27410829e3
    compute_wf = ".false."


@ex.capture
def get_id(_run):
    return _run._id


@ex.automain
def run(d, n, n_remove, Z1, Z2, m, c, R, ximax, ximin, epsilon, eta_slp, xi_slp, save_step, tot_diag, maxit, eig, compute_wf):
    # Create temporary folder for result
    result_folder = f"tmp_{get_id()}"
    try:
        os.mkdir(result_folder)
    except Exception as error:
        print(f"Failed to create output directory: {error}")
        raise
    print(f"Temporary folder created: {result_folder}")

    # Write formatted input file
    with open(f"{result_folder}/input.txt", "w") as f:
        f.write(f"""\
    {d},                   - d : Order of the B-Spline (order Mathematica + 1)
    {n},                  - n : Number of used B-Splines
    {n_remove},                   - n_remove : Number of B-Splines to remove at each end
    {Z1:.2f},                - Z1 : Z1 value
    {Z2:.2f},                - Z2 : Z2 value
    {m:.2f},                - m : electron mass
    {c:.15e},       - c : speed of light
    {R:.15e},                 - R : SEMI-interatomic distance
    {ximax:.15e},                - ximax : Maximum xi value
    {ximin:.15e},                 - ximin : Minimum xi value
    {epsilon:.15e},                 - epsilon : to avoid singularities
    {eta_slp:.15e},              - eta_slp : slope value to generate knots on eta
    {xi_slp:.15e},              - xi_slp : slope value to generate knots on xi
    {save_step},              - save_step : Save every matrices
    {result_folder},              - name of the temporary folder
    {tot_diag},              - tot_diag : Perform total diagonalization (true) or partial (false)
    {maxit},              - maxit : Maximum number of iterations for the eigensolver
    {eig:.15e},              - eig : Initial eigenvalue guess
    {compute_wf},              - compute_wf : Compute wavefunction (true/false)
""")

    # Compile and run Fortran program
    print("Compiling and running Fortran code...")
    subprocess.run(["make", "-C", "fortran", "main.out"],
                   check=True, text=True)
    # result = subprocess.run(["./fortran/main.out"], input=open("input.txt").read(), text=True, capture_output=True)
    result = subprocess.Popen(["./fortran/main.out"], stdin=open(f"{result_folder}/input.txt"),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Monitor memory and CPU usage
    pid = result.pid
    p = psutil.Process(pid)
    print(f"Monitoring process PID: {pid}")
    # p.cpu_percent(interval=None) # Initialize CPU percent calculation
    # try:
    #     while True:
    #         ex.log_scalar("memory_mb", p.memory_info().rss / (1024 * 1024))
    #         ex.log_scalar("cpu_percent", p.cpu_percent(interval=None))
    #         if result.poll() is not None:
    #             break
    #         time.sleep(5)
    # except psutil.NoSuchProcess:
    #     pass

    # Save raw logs and outputs
    with open(f"{result_folder}/output.log", "w") as f:
        f.write(result.stdout.read())

    ex.add_artifact(f"{result_folder}/input.txt")
    ex.add_artifact(f"{result_folder}/output.log")
    ex.add_artifact(f"{result_folder}/eigenvalues.txt")
    ex.add_artifact(f"{result_folder}/log_file")

    if compute_wf.lower() == ".true.":
        try:
            ex.add_artifact(f"{result_folder}/wavefun.txt")
        except Exception as error:
            print(f"Failed to add wavefun.txt artifact: {error}")

    for file in glob.glob(f"{result_folder}/*.csv"):
        ex.add_artifact(file)

    branch = subprocess.check_output(
        ["git", "-C", "fortran", "rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
    cwd = subprocess.check_output(["pwd"]).decode().strip()
    msg = f"✅ Computation finished on branch \"{branch}\" in directory \"{cwd}\""
    push(msg)

    last_eigenvalue = extract_last_eigenvalue(
        f"{result_folder}/eigenvalues.txt")

    # try:
    #     os.remove(f"{result_folder}/input.txt")
    #     os.remove(f"{result_folder}/output.log")
    #     os.remove(f"{result_folder}/eigenvalues.txt")
    #     os.remove(f"{result_folder}/log_file")
    # except Exception as error:
    #     print(error)
    #
    # for file in glob.glob(f"{result_folder}/*.csv"):
    #     try:
    #         os.remove(file)
    #     except Exception as error:
    #         print(error)

    return {"last_eigenvalue": last_eigenvalue}
