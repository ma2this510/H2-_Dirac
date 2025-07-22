from sacred import Experiment
from sacred.observers import FileStorageObserver, MongoObserver

import os

# Setup Sacred experiment
ex = Experiment("manual_entry")
ex.observers.append(FileStorageObserver("experiments"))  # logs stored in ./experiments/
ex.observers.append(MongoObserver.create())

# Input parameters
@ex.config
def cfg():
    d = 8
    n = 20
    n_remove = 0
    Z1 = 1.00
    Z2 = 1.00
    m = 1.00
    c = 137.035999679
    R = 1.0
    ximax = 30.0
    ximin = 1.0
    epsilon = 0.0
    eta_slp = 4.0e-2
    save_step = ".false."

    # Required output files
    eigen_file = "eigenvalues.txt"
    C11_file = "C11one.csv"
    C22_file = "C22one.csv"

@ex.automain
def record_manual_run(eigen_file, C11_file, C22_file):
    # Ensure files exist and register them as artifacts
    for f in [eigen_file, C11_file, C22_file]:
        if os.path.exists(f):
            ex.add_artifact(f)
        else:
            print(f"⚠️ Warning: File not found — {f}")

    print("✅ Manual result entry completed.")

