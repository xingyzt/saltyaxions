import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import LinearNDInterpolator
import multiprocessing

plt.rcParams['figure.figsize'] = [16, 9]
plt.rcParams['figure.dpi'] = 300


# plt.rcParams.update({
#     "lines.color": "white",
#     "patch.edgecolor": "white",
#     "text.color": "black",
#     "axes.facecolor": "white",
#     "axes.edgecolor": "lightgray",
#     "axes.labelcolor": "white",
#     "xtick.color": "white",
#     "ytick.color": "white",
#     "grid.color": "lightgray",
#     "figure.facecolor": "black",
#     "figure.edgecolor": "black",
#     "savefig.facecolor": "black",
#     "savefig.edgecolor": "black"
# })


isotopes = [
    'c12', 
    'o16', 
    'ne20', 
    'na23',
]


names = {
    'c12': '$^{12}$C', 
    'o16': '$^{16}$O', 
    'ne20': '$^{20}$Ne', 
    'na23': '$^{23}$Na',
    }

color = {
    'c12': 'skyblue',
    'o16': 'violet',
    'ne20': 'springgreen',
    'na23': 'coral',
}

widths = {
    'c12': 2,
    'o16': 2,
    'ne20': 2,
    'na23': 4,
}

m = 9
g = -10
key = f'm{m:04.1f}_g{g:+04.2f}_0'
index = pd.read_csv(f'csv/{key}/index.csv')
t = np.array(index["time to core O depletion (years)"])
N = np.array(index["surface num_a (/s)"]) * 100
slices = [ pd.read_csv(f'csv/{key}/slice_{int(i)}.csv') for i in index["slice"]]

fig, (tax, Xax) = plt.subplots(nrows=2, figsize=(0.6*16, 0.6*9), dpi=300, gridspec_kw={'height_ratios': [2, 5]})
Tax = Xax.twinx()

for (j, jt) in enumerate(list(10**np.linspace(7, -1, 30*60))):

    i = np.argwhere(t < jt)[0][0]
    s = slices[i]

    plt.sca(tax)
    tax.set_yscale('log')
    tax.set_xscale('log')
    plt.plot(t, N, c='black')
    plt.axvline(t[i], c='blue')
    plt.ylim(1e43, 1e47)
    plt.xlim(1e7, 1e-1)
    plt.title("Axion flux from a 9 $M_\odot$ star with $g_{eff} = 10^{-10}$", c="white")
    plt.xlabel("Years to core collapse")
    plt.ylabel("Axions / sec")
    plt.grid(c="gray", ls=":")
    
    plt.sca(Xax)
    Xax.set_yscale('log')

    m = np.array(s["mass enclosed (Msun)"])

    for iso in isotopes:
        plt.plot(m, np.array(s["X_"+iso]), label=names[iso], c=color[iso], lw=widths[iso])

    plt.ylim(1e-4, 1)
    plt.xlim(0, 9)
    plt.xlabel("Mass enclosed ($M_\odot$)")
    plt.ylabel("Isotope mass fraction")
    plt.grid(c="gray", ls=":")
    plt.legend()

    plt.sca(Tax)
    Tax.set_yscale('log')
    plt.plot(m, np.array(s["T (K)"]), ls="--", c='black')
    plt.ylim(1e6, 1e10)
    plt.ylabel("Temperature (K)")
    plt.tight_layout()


    plt.savefig(f'slice2/{j:05}.png')
    tax.clear()
    Xax.clear()
    Tax.clear()
    