from sacred import Experiment
from sacred.observers import FileStorageObserver, MongoObserver
from dotenv import load_dotenv
import subprocess
import os
import glob

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
    subprocess.run(["make", "-C", "fortran", "main.out"], check=True)
    result = subprocess.run(["./fortran/main.out"], input=open("input.txt").read(), text=True, capture_output=True)

    # Save raw logs and outputs
    with open("output.log", "w") as f:
        f.write(result.stdout)

    ex.add_artifact("input.txt")
    ex.add_artifact("output.log")
    ex.add_artifact("eigenvalues.txt")
    ex.add_artifact("log_file")

    for file in glob.glob("*.csv"):
        ex.add_artifact(file)
