import numpy as np
from concurrent.futures import ThreadPoolExecutor
import subprocess
import threading
import time
import sys

thread_num = int(sys.argv[1]) if len(sys.argv) > 1 else 1 

worker_counter = 0
counter_lock = threading.Lock()

xi_slp_list = np.round(np.linspace(10,11, 10), 7)
eta_slp_list = np.round(np.linspace(0.93, 0.97, 12), 7)

param_list = np.array(np.meshgrid(eta_slp_list, xi_slp_list)).T.reshape(-1, 2)

def run_command(params):
    eta_slp, xi_slp = params
    command = f"python3 run_experiment.py with n=26 d=10 ximax=30 eta_slp={eta_slp} xi_slp={xi_slp} -c 'analysis grid slp param'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    print("Subprocess finished with return code:", result.returncode)
    if result.returncode != 0:
        print("Subprocess error output:", result.stderr)
        return np.inf
    output = result.stdout
    
    value = np.inf
    for line in output.splitlines():
        if "Last eigenvalue extracted:" in line:
            # Extract the float part
            num = line.split(":")[1].strip()
            value = np.float64(num)
            break

    print(f"Evaluated parameters: xi_slp={xi_slp}, eta_slp={eta_slp}, xi_max={xi_max} => last_eigenvalue={value}")
    
    e_ref = -1.10264158103257716412
    error = abs(value - e_ref)
    print(f"Error with respect to reference: {error}")
    return error

def delayed_run_fun(params):
    global worker_counter
    with counter_lock:
        worker_id = worker_counter
        worker_counter += 1

    delay = (worker_id % thread_num) * 1.0  # 1 sec delay between workers
    time.sleep(delay)
    print(f"Worker {worker_id} starting after {delay:.1f}s delay")
    return run_command(params)

with ThreadPoolExecutor(max_workers=thread_num) as executor:
    executor.map(delayed_run_fun, param_list)
