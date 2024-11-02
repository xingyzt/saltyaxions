import csv
import os
import subprocess

if not os.path.exists('sims'):
    os.makedirs('sims')

if not os.path.exists('logs'):
    os.makedirs('logs')

m = float(input("mass: m/Msun = "))
g = float(input("coupling: log(g) = "))

template_path = 'mesa'
path = f'm{m:04.1f}_g{g:+04.2f}'
sim_path = f'sims/{path}'
log_path = f'logs/{path}'

subprocess.run(['cp', '-r', template_path, sim_path])

with open(f'{template_path}/inlist_make_late_pre_zams') as fin:
    with open(f'{sim_path}/inlist_make_late_pre_zams', "w") as fout:
        fout.write(
            fin.read().replace(f'initial_mass = 20.0', f'initial_mass = {m:.9f}')
        )

with open(f'{template_path}/inlist_common') as fin:
    with open(f'{sim_path}/inlist_common', "w") as fout:
        fout.write(
            fin.read().replace(f'x_ctrl(99) = 1e-09', f'x_ctrl(99) = {10**g:.9e}')
        )

os.system(f'cd {sim_path} && ./rn_all && cd ..')
os.system(f'mv {sim_path}/LOGS {log_path}')
