import mesa_reader as mr
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import csv
import os

slices = False

masses = [
    15, 16, 17, 18, 19, 20, 21, 22
]

couplings = [
    -8.5
]

isotopes = [
    'c12', 
    'o16', 
    'ne20', 
    'na23',
]

labels = {
    'model': 'slice',
    'age': 'age (years)',
    'til': 'time to core O depletion (years)',
    'dt': 'dt (years)',
    
    'm': 'mass enclosed (Msun)',
    'dm': 'dm (g)',
    'r': 'r (Rsun)',
    'dr': 'dr (cm)',
    'T': 'T (K)',
    'T_core': 'core T (K)',
    'T_eff': 'effective T (K)',
    'log_T': 'log T',
    
    'eps_grav': 'eps_grav (ergs/g s)',
    'eps_nuc': 'eps_nuc (ergs/g s)', 
    'eps_non_nuc_neu': 'eps_neu (ergs/g s)',
    'eps_a': 'eps_a (ergs/g s)',
    
    'lum_gamma': 'lum_gamma (ergs/s)',
    'lum_neu': 'lum_neu (ergs/s)',
    'lum_a': 'lum_a (ergs/s)',
    
    'lum_gamma_surf': 'surface lum_gamma (ergs/s)',
    'lum_neu_surf': 'surface lum_neu (ergs/s)',
    'lum_a_surf': 'surface lum_a (ergs/s)',
}

for iso in isotopes:
    labels[iso] = 'X_' + iso
    labels['log_' + iso] = 'log X_' + iso
    labels['X_' + iso] = 'avg X_' + iso

def get_profile(i):

    l = mr.MesaLogDir(inpath)
    #bulk
    b = mr.MesaData(f'{inpath}/profile{l.profile_numbers[i]}.data')

    #header
    h = b.header_data 

    # print(h)

    # profile
    p = {
        'profile': l.profile_numbers[i],
        'model': l.model_numbers[i],
        'age': h['star_age'],
        'T_eff': h['Teff'],
        
        'T': b.data('T'), # Kelvins,
        'log_T': b.data('logT'),
        'dm': b.data('dm'), # g
        'm': b.data('mass'), # M_sun
        'dt': h['time_step'], # year
        'r': b.data('radius'), # R_sun
        'dr': b.data('dr'), # cm
        'eps_grav': b.data('eps_grav'), # ergs / g s
        'eps_nuc': b.data('eps_nuc'), 
        # 'eps_nuc_minus_non_nuc_neu': b.data('eps_nuc_minus_non_nuc_neu'),
        # 'eps_nuc_neu': b.data('eps_nuc_neu_total'),
        # 'eps_nuc_plus_nuc_neu': b.data('eps_nuc_plus_nuc_neu'),
        # 'eps_nuc_start': b.data('eps_nuc_start'),
        'eps_non_nuc_neu': b.data('non_nuc_neu'),
        'eps_a': b.data('axion'),
        # 'net_nuclear_energy': 10**b.data('net_nuclear_energy'),
        # 'net_energy': 10**b.data('net_energy'),
        'lum': b.data('lum_erg_s')
    }

    for iso in isotopes:
        p[ iso ] = b.data( iso )
        p[ 'log_' + iso ] = np.nan_to_num(b.data( 'log_' + iso ), -99)

        p[ 'X_' + iso ] = np.sum(p['dm'] * p[iso]) / np.sum(p['dm'])

    p['lum_a'] = np.cumsum((p['eps_a'] * p['dm'])[::-1])[::-1]
    p['lum_neu'] = np.cumsum((p['eps_non_nuc_neu'] * p['dm'])[::-1])[::-1]
    p['lum_gamma'] = np.cumsum((p['eps_nuc'] * p['dm'])[::-1])[::-1]

    p['lum_a_surf'] = p['lum_a'][-1]
    p['lum_neu_surf'] = p['lum_neu'][-1]
    p['lum_gamma_surf'] = p['lum_gamma'][-1]
    p['T_core'] = p['T'][0]

    return p

for m in masses: 
    for g in couplings:

        path = f'm{m:.1f}_g{g:+.2f}'
        inpath = f'mesa/{path}'
        outpath = f'csv/{path}'
        print(path)

        l = mr.MesaLogDir(inpath)
        N = len(l.profile_numbers)
        # N = 2
        profiles = list(ProcessPoolExecutor().map(get_profile, range(N)))


        if not os.path.exists(outpath):
            os.makedirs(outpath)

        for p in profiles:
            p['til'] = profiles[N-1]['age'] - p['age']

        with open(f'{outpath}/index.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            p = profiles[0]
            writer.writerow([ 
                labels[key] 
                for key in labels 
                if np.isscalar(p[key])
            ])
            for p in profiles:
                writer.writerow(np.array([ 
                    p[key]
                    for key in labels
                    if np.isscalar(p[key])
                ]).T)

        if not slices: 
            continue

        for p in profiles:

            with open(f'{outpath}/slice_{p["model"]}.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([ 
                    labels[key] 
                    for key in labels 
                    if not np.isscalar(p[key])
                ])
                writer.writerows(np.array([ 
                    p[key][::-1]
                    for key in labels
                    if not np.isscalar(p[key])
                ]).T)
