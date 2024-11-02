# Stellar axion simulation

## Simulating with MESA

- We use MESA version 24.03.1. You can follow the official instructions to install MESA.
- Then, clone this repository and navigate to the `mesa` directory. Run the `clean` and `mk` scripts to build the `star` executable.
- Then, run `simulate.py`.

## Data processing

- To convert your MESA logs to csv format, 
  install numpy, pandas, and [py_mesa_reader](https://github.com/wmwolf/py_mesa_reader) for Python,
  then run this repository's `reduce.py`.

