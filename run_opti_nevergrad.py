import nevergrad as ng
from concurrent import futures
import numpy as np
import subprocess
import re
import time
import threading

worker_counter = 0
counter_lock = threading.Lock()

thread_num = 24 
max_run = 700

print("Nevergrad version:", ng.__version__)
print("Numpy version:", np.__version__)

def run_fun(xi_slp, eta_slp, xi_max):

    print(f"Running with parameters: xi_slp={xi_slp}, eta_slp={eta_slp}, xi_max={xi_max}")

    command = f"python3 run_experiment.py with n=26 d=10 ximax={xi_max} eta_slp={eta_slp} xi_slp={xi_slp} -c 'Nevergrad optimization test 6'"

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

def delayed_run_fun(xi_slp, eta_slp, xi_max):
    global worker_counter
    with counter_lock:
        worker_id = worker_counter
        worker_counter += 1

    delay = (worker_id % thread_num) * 2.0  # 1 sec delay between workers
    time.sleep(delay)
    print(f"Worker {worker_id} starting after {delay:.1f}s delay")
    return run_fun(xi_slp, eta_slp, xi_max)

instrum = ng.p.Instrumentation(
    ng.p.Scalar(lower=0, upper=10),  # xi_slp
    ng.p.Scalar(lower=0, upper=1),  # eta_slp
    ng.p.Scalar(lower=1, upper=100)  # xi_max
)

optimizer = ng.optimizers.registry["TwoPointsDE"](parametrization=instrum, budget=max_run, num_workers=thread_num)
with futures.ThreadPoolExecutor(max_workers=optimizer.num_workers) as executor:
    recommendation = optimizer.minimize(delayed_run_fun, executor=executor, batch_mode=False)

print("Best parameters found: ", recommendation.value)
