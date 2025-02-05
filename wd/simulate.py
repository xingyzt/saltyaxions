import csv
import os
import subprocess

if not os.path.exists('sims'):
    os.makedirs('sims')

if not os.path.exists('logs'):
    os.makedirs('logs')

m = float(input("mass: m/Msun = "))
log_g = float(input("coupling: log(g) = "))
affect = bool(int(input("evolution weight (1 or 0): ")))
g = 10**log_g

template_path = 'mesa'
path = f'm{m:04.1f}_g{log_g:+05.2f}{"" if affect else "_0"}'
sim_path = f'sims/{path}'
log_path = f'logs/{path}'

subprocess.run(['cp', '-r', template_path, sim_path])

with open(f'{template_path}/inlist_common') as fin:
    with open(f'{sim_path}/inlist_common', "w") as fout:
        fout.write(
            fin.read().replace(f'initial_mass = 11.0', f'initial_mass = {m:.9f}')
            .replace(f'x_ctrl(99) = 1e-09', f'x_ctrl(99) = {g:.9e}')
            .replace(f'x_logical_ctrl(99) = .true.', f'x_logical_ctrl(99) = {".true" if affect else ".false."}')
        )

os.system(f'cd {sim_path} && ./rn_all && cd ..')
os.system(f'mv {sim_path}/LOGS {log_path}')
