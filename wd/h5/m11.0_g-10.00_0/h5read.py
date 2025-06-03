import h5py
import matplotlib.pyplot as plt
with h5py.File('index.hdf5', 'r') as f:
    indices = f['slice'][()]
with h5py.File('slices.hdf5', 'r') as f:
    print(f[f'slice_1'].keys())
    for i in indices:
        x = f[f'slice_{i}']['r (Rsun)'][:]
        y = f[f'slice_{i}']['T (K)'][:]
        plt.plot(x, y)
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.xlabel('r (Rsun)')
    plt.ylabel('T (K)')
    plt.show()