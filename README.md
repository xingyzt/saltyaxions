# Stellar axion simulation

## Simulating with MESA

- We use MESA version 24.03.1. You can follow the official instructions to install MESA.
- Then, clone this repository and navigate to the `mesa` directory. Run the `clean` and `mk` scripts to build the `star` executable.
- To change the ZAMS mass, edit the `initial_mass` parameter in the `&controls` section of `mesa/inlist_make_late_pre_zams`.
  - If you changed the ZAMS mass, then you should start simulating by running the `rn_all` script. 
  If you didn't, then you can run with the `rn_from_c_burn` script if MESA saved a `after_core_he_burn.mod` from a previous run. 
  The axion losses due to Na23 are insignificant before the start of core C burning.
- To change the effective axion-nuclear coupling, edit the `x_ctrl(99)` parameter in the `&controls` section of `mesa/inlist_common`.
  - The axion losses are calculated in `mesa/src/run_star_extras.f90`.
- After a simulation with `[M]` solar masses and coupling `10^[G]`, you should rename the `mesa/LOGS` directory to `logs/m[M]_g[G]`,
  where `[M]` is unsigned with one decimal place, and `[G]` is signed with two decimal places. For example, `logs/m20.0_g-8.50`.

## Data processing

- To convert your MESA logs to csv format, 
  install numpy, pandas, and [py_mesa_reader](https://github.com/wmwolf/py_mesa_reader) for Python,
  then run this repository's `reduce.py`.

