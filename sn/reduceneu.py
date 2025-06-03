import mesa_reader as mr
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import csv
import os

slices = False
affect = False
masses = [ 10 ]
couplings = [ -10 ]


isotopes = [
    'c12', 
    'o16', 
    'ne20', 
    'na23',
]

neu_types = [ 'nuc', 'non_nuc' ]
non_nuc_neu_subtypes = [ 'plas', 'brem', 'phot', 'pair', 'reco' ]
gamma_types = [ 'nuc', 'grav' ]

labels = {
    'model': 'slice',
    'm': 'mass (Msun)',
    'coupling': 'coupling',
    'age': 'age (years)',
    'til': 'time to core O depletion (years)',
    'dt': 'dt (years)',
    
    'm_enc': 'mass enclosed (Msun)',
    'dm': 'dm (g)',
    'r': 'r (Rsun)',
    'dr': 'dr (cm)',
    'T': 'T (K)',
    'T_core': 'core T (K)',
    'T_eff': 'effective T (K)',
    'log_T': 'log T',
    
    'eps_gamma': 'eps_gamma (ergs/g s)', 
    'eps_neu': 'eps_neu (ergs/g s)',
    'eps_a': 'eps_a (ergs/g s)',
    
    'lum_gamma': 'lum_gamma (ergs/s)',
    'lum_neu': 'lum_neu (ergs/s)',
    'num_a': 'num_a (/s)',
    
    'lum_gamma_surf': 'surface lum_gamma (ergs/s)',
    'lum_neu_surf': 'surface lum_neu (ergs/s)',
    'lum_a_surf': 'surface lum_a (ergs/s)',
    'num_a_surf': 'surface num_a (/s)',
    
    'cum_e_gamma': 'cumulative e_gamma (ergs)',
    'cum_e_neu': 'cumulative e_neu (ergs)',
    'cum_e_a': 'cumulative e_a (ergs)',
}

for iso in isotopes:
    labels['X_' + iso] = 'X_' + iso
    labels['log_X_' + iso] = 'log X_' + iso
    labels['avg_X_' + iso] = 'avg X_' + iso

for neu_type in [*neu_types, *non_nuc_neu_subtypes]:
    labels[f'lum_{neu_type}_neu_surf'] = f'surface lum_{neu_type}_neu (ergs/s)'

for gamma_type in gamma_types:
    labels[f'lum_{gamma_type}_gamma_surf'] = f'surface lum_{gamma_type}_gamma (ergs/s)'

labels['lum_mesa_gamma_surf'] = 'MESA surface lum_gamma (ergs/s)'

def eps_a(na23, T, g):
     eps0 = 8.6E+27 # ergs / g s; baseline
     T0 = 5.106E+9 # Kelvins; baseline
     mu0 = 1.5 # unitless; chemical potential

     return eps0 * na23 * g * g * np.exp(-T0/T) / ( 1 + mu0 * np.exp(-T0/T) )

def get_profile(i):

    l = mr.MesaLogDir(inpath)
    #bulk
    b = mr.MesaData(f'{inpath}/profile{l.profile_numbers[i]}.data')

    #header
    h = b.header_data 

    # profile
    p = {
        'profile': l.profile_numbers[i],
        'model': l.model_numbers[i],
        'age': h['star_age'],
        'T_eff': h['Teff'],
        'g_eff': h['g_eff'],
        'lum_a_surf': h['lum_axion_surf'],
        
        'T': b.data('T'), # Kelvins,
        'log_T': b.data('logT'),
        'dm': b.data('dm'), # g
        'm_enc': b.data('mass'), # M_sun
        'dt': h['time_step'], # year
        'r': b.data('radius'), # R_sun
        'dr': b.data('dr'), # cm
        'eps_nuc': b.data('eps_nuc'), 
        'eps_nuc_minus_non_nuc_neu': b.data('eps_nuc_minus_non_nuc_neu'),
        # 'eps_nuc_start': b.data('eps_nuc_start'),
        # 'eps_nuc_plus_nuc_neu': b.data('eps_nuc_plus_nuc_neu'),
        'eps_nuc_neu': b.data('eps_nuc_neu_total'),
        'eps_non_nuc_neu': b.data('non_nuc_neu'),
        'eps_grav_gamma': b.data('eps_grav'), # ergs / g s
        'eps_nuc_gamma': b.data('eps_nuc'), # ergs / g s
        'eps_gamma': b.data('net_energy'),
        'eps_a': b.data('axion'),
        # 'net_nuclear_energy': 10**b.data('net_nuclear_energy'),
        # 'net_energy': 10**b.data('net_energy'),
        'lum_mesa_gamma': b.data('lum_erg_s')
    }

    # data slices: index 0 is surface, index -1 is center

    for iso in isotopes:
        p[ 'X_' + iso ] = b.data( iso )
        p[ 'log_X_' + iso ] = np.nan_to_num(b.data( 'log_' + iso ), -99)

        p[ 'avg_X_' + iso ] = np.sum(p['dm'] * p[ 'X_' + iso]) / np.sum(p['dm'])

    p['eps_a'] = eps_a(na23=p['X_na23'], T=p['T'], g=p['g_eff'])
    p['lum_a'] = np.cumsum((p['eps_a'] * p['dm'])[::-1])[::-1]
    p['lum_a_surf'] = p['lum_a'][0]
    p['num_a'] = p['lum_a'] * 1418524 # axions/erg
    p['num_a_surf'] = p['num_a'][0]

    for neu_type in neu_types:
        p[f'lum_{neu_type}_neu'] = np.cumsum((p[f'eps_{neu_type}_neu'] * p['dm'])[::-1])[::-1]
        p[f'lum_{neu_type}_neu_surf'] = p[f'lum_{neu_type}_neu'][0] 

    for neu_type in non_nuc_neu_subtypes:
        p[f'eps_{neu_type}_neu'] = b.data(f'nonnucneu_{neu_type}')
        p[f'lum_{neu_type}_neu'] = np.cumsum((p[f'eps_{neu_type}_neu'] * p['dm'])[::-1])[::-1]
        p[f'lum_{neu_type}_neu_surf'] = p[f'lum_{neu_type}_neu'][0]

    p['eps_neu'] = p['eps_non_nuc_neu'] + p['eps_nuc_neu']
    p['lum_neu'] = p['lum_non_nuc_neu'] + p['lum_nuc_neu']
    p['lum_neu_surf'] = p['lum_neu'][0]

    for gamma_type in gamma_types:
        p[f'lum_{gamma_type}_gamma'] = np.cumsum((p[f'eps_{gamma_type}_gamma'] * p['dm'])[::-1])[::-1]
        p[f'lum_{gamma_type}_gamma_surf'] = p[f'lum_{gamma_type}_gamma'][0] 

    p['lum_gamma'] = np.cumsum((p['eps_gamma'] * p['dm'])[::-1])[::-1]
    p['lum_gamma_surf'] = p['lum_gamma'][0]

    p['lum_mesa_gamma_surf'] = p['lum_mesa_gamma'][0]

    p['T_core'] = p['T'][-1]

    p['m'] = p['m_enc'][0]

    return p

for m in masses: 
    for g in couplings:

        path = f'm{m:04.1f}_g{g:+05.2f}{"" if affect else "_0"}'
        inpath = f'logs/{path}'
        outpath = f'csv/{path}'
        print(path)

        l = mr.MesaLogDir(inpath)
        N_raw = len(l.profile_numbers)
        profiles_raw = list(ProcessPoolExecutor().map(get_profile, range(N_raw)))
        profiles = []

        dt = 0
        for (i, p) in enumerate(profiles_raw):
            dt += p['dt']
            if dt > 1e-6: # threshold for minimum timestep
                p['dt'] = dt
                profiles.append(p)
                dt = 0
        N = len(profiles)

        if not os.path.exists(outpath):
            os.makedirs(outpath)

        for (i, p) in enumerate(profiles):

            p['til'] = profiles[N-1]['age'] - p['age']

            p['coupling'] = g

            p['cum_e_a'] = p['lum_a_surf'] * p['dt'] * 31536000
            p['cum_e_neu'] = p['lum_neu_surf'] * p['dt'] * 31536000
            p['cum_e_gamma'] = p['lum_gamma_surf'] * p['dt'] * 31536000
            
            if i > 0:
                p['cum_e_a'] += profiles[i-1]['cum_e_a']
                p['cum_e_neu'] += profiles[i-1]['cum_e_neu']
                p['cum_e_gamma'] += profiles[i-1]['cum_e_gamma']

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
