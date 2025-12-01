import optuna
import subprocess
import numpy as np
import threading
import time
# from pymongo import MongoClient
# from dotenv import load_dotenv
# import pandas as pd
# import os

# load_dotenv()

thread_num = 8 
trial_num = 300
n = 12

comment_id = f"Optuna optimization test n={n} default"

worker_counter = 0
counter_lock = threading.Lock()

def run_fun_optuna(trial):
    xi_slp_norm = trial.suggest_float("xi_slp_norm", 0.0, 1.0)
    eta_slp_norm = trial.suggest_float("eta_slp_norm", 0.0, 1.0)
    xi_max_norm = trial.suggest_float("xi_max_norm", 0.0, 1.0)

    # Wait to prevent main.out being overloaded
    global worker_counter
    with counter_lock:
        worker_id = worker_counter
        worker_counter += 1
    delay = (worker_id % thread_num) * 1.2  # 1.2 sec delay between workers
    print(f"Worker {worker_id} starting after {delay:.1f}s delay")
    time.sleep(delay)

    # Scale parameters to their actual ranges
    xi_slp = xi_slp_norm * 10.0                  # 0 to 10
    eta_slp = eta_slp_norm * 1.0              # 0 to 1 
    xi_max = xi_max_norm * (100.0 - 1.0) + 1.0  # 1 to 100

    print(f"Running with parameters: xi_slp={xi_slp}, eta_slp={eta_slp}, xi_max={xi_max}")

    command = f"python3 run_experiment.py with n={n} d=10 ximax={xi_max} eta_slp={eta_slp} xi_slp={xi_slp} -c '{comment_id}'"

    result = subprocess.run(command, shell=True, capture_output=True)
    print("Subprocess finished with return code:", result.returncode)
    if result.returncode != 0:
        print("Subprocess error output:", result.stderr)
        return np.inf
    output = result.stdout
    
    value = np.inf
    for line in output.splitlines():
        try:
            if b"Last eigenvalue extracted:" in line:
                # Extract the float part
                num = str(line).split(":")[1].replace("'", "").strip()
                value = np.float128(num)
        except Exception as err:
            print(f"Unexpected error : {err}")
            break

    print(f"Evaluated parameters: xi_slp={xi_slp}, eta_slp={eta_slp}, xi_max={xi_max} => last_eigenvalue={value}")
    
    e_ref = np.float128(-1.10264158103257716412)
    error = np.abs(value - e_ref)
    print(f"Error with respect to reference: {error}")
    return np.log10(error)


study = optuna.create_study()
study.optimize(run_fun_optuna, n_trials=trial_num, n_jobs=thread_num)

best_params = study.best_params
print("Best parameters found (norm):", best_params)
# Scale back to actual ranges
best_xi_slp = best_params["xi_slp_norm"] * 10.0
best_eta_slp = best_params["eta_slp_norm"] * 1.0
best_xi_max = best_params["xi_max_norm"] * (100.0 - 1.0) + 1.0
print(f"Best parameters found (actual): xi_slp={best_xi_slp}, eta_slp={best_eta_slp}, xi_max={best_xi_max}")