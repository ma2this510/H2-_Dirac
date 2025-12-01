import nevergrad as ng
from concurrent import futures
import numpy as np
import subprocess
import re
import time
import threading
from pymongo import MongoClient
from dotenv import load_dotenv
import pandas as pd
import os

worker_counter = 0
counter_lock = threading.Lock()

thread_num = 8 
max_run = 300

n = 12

comment_id = f"Nevergrad optimization test n={n} CMA"

print("Nevergrad version:", ng.__version__)
print("Numpy version:", np.__version__)

def run_fun(xi_slp_norm, eta_slp_norm, xi_max_norm):

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

def delayed_run_fun(xi_slp_norm, eta_slp_norm, xi_max_norm):
    global worker_counter
    with counter_lock:
        worker_id = worker_counter
        worker_counter += 1

    delay = (worker_id % thread_num) * 2.0  # 1 sec delay between workers
    time.sleep(delay)
    print(f"Worker {worker_id} starting after {delay:.1f}s delay")
    return run_fun(xi_slp_norm, eta_slp_norm, xi_max_norm)

instrum = ng.p.Instrumentation(
    ng.p.Scalar(lower=0, upper=1),  # xi_slp_norm
    ng.p.Scalar(lower=0, upper=1),  # eta_slp_norm
    ng.p.Scalar(lower=0, upper=1)  # xi_max_norm
)

optimizer = ng.optimizers.CMA(parametrization=instrum, budget=max_run, num_workers=thread_num)

# -------------------------------------------------------------------------------------------------
# Import already existing data from MongoDB
load_dotenv()

client = MongoClient(os.environ['URI'])
db = client[os.environ['DB']]
runs_collection = db['runs']

runs = list(runs_collection.find({"meta.comment": comment_id}, {"config": 1, "result": 1}))
if not runs:
    print("No completed runs found in the given range.")
    
else:
    df = pd.DataFrame(runs)
    print(df.count())

    df_conf = pd.json_normalize(df['config'])
    df_conf.columns = [f'conf_{col}' for col in df_conf.columns]

    df_res = pd.json_normalize(df['result'])
    df_res.columns = [f'res_{col}' for col in df_res.columns]

    col = df.columns.difference(['config','result'])
    df_final_01 = pd.concat([df[col], df_conf, df_res], axis=1)

    df_final_01['res_last_eigenvalue'] = df_final_01['res_last_eigenvalue'].astype(np.float64)

    e_ref = -1.10264158103257716412

    df_final_01 = df_final_01.sort_values(by=['_id']).reset_index(drop=True)

    df_final_01['res_error'] = abs(df_final_01['res_last_eigenvalue'] - e_ref)
    df_final_01['res_log_error'] = np.log10(df_final_01['res_error'])

    for index, row in df_final_01.iterrows():
        xi_slp = row['conf_xi_slp']
        eta_slp = row['conf_eta_slp']
        xi_max = row['conf_ximax']
        log_error = row['res_log_error']

        xi_slp_norm = xi_slp / 10.0
        eta_slp_norm = eta_slp / 1.0
        xi_max_norm = (xi_max - 1.0) / (100.0 - 1.0)

        #candidate = optimizer.parametrization.spawn_child(new_value={xi_slp, eta_slp, xi_max})
        optimizer.suggest(xi_slp_norm, eta_slp_norm, xi_max_norm)
        candidate = optimizer.ask()
        optimizer.tell(candidate, log_error)
        print(f"Imported run {index}: xi_slp={xi_slp}, eta_slp={eta_slp}, xi_max={xi_max}, log_error={log_error}")

# -------------------------------------------------------------------------------------------------

with futures.ThreadPoolExecutor(max_workers=optimizer.num_workers) as executor:
    recommendation = optimizer.minimize(delayed_run_fun, executor=executor, batch_mode=True)

print("Best parameters found: ", recommendation.value)
