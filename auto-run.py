import numpy as np
import subprocess
import time
import psutil
import os

mem = psutil.virtual_memory().percent
cpu = psutil.cpu_percent(interval=1)

# eta_slp_list = np.round(np.linspace(5, 9, 10), 2)
# xi_slp_list = np.round(np.linspace(0.6, 0.95, 10), 2)
eta_slp_list = 0.82
xi_slp_list = np.round(np.linspace(7.6, 8.6, 40), 2)
param_list = np.array(np.meshgrid(eta_slp_list, xi_slp_list)).T.reshape(-1, 2)

if os.path.exists("current_run"):
    try:
        with open("current_run", "r") as f:
            i_run = int(f.read().strip())
    except Exception:
        i_run = 0
        with open("current_run", "w") as f:
            f.write("0")
else:
    i_run = 0
    with open("current_run", "w") as f:
        f.write("0")

with open("current_run", "w") as f:
    while i_run < param_list.shape[0]:
        mem = psutil.virtual_memory().percent
        cpu = psutil.cpu_percent(interval=None)

        if mem < 50 and cpu < 50 :
            eta_slp = param_list[i_run, 0]
            xi_slp = param_list[i_run, 1]

            command = f"nohup python3 run_experiment.py with n=24 d=10 ximax=30 eta_slp={eta_slp} xi_slp={xi_slp} -c 'test xi regime quarter 3' &"

            subprocess.run(command, shell=True)
            i_run += 1
            f.seek(0)
            f.write(str(i_run))
            print(f"Submitted run {i_run}/{param_list.shape[0]} with eta_slp={eta_slp} and xi_slp={xi_slp}")
            
            if i_run % 5 == 0:
                time.sleep(120)
        else:
            print(f"System busy (Memory: {mem}%, CPU: {cpu}%), waiting to submit next job...")
            time.sleep(400)

        time.sleep(20)
