#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:37:42 2022

@author: jtm545

Script to make the stimuli for the MonBin1 experiment. 

"""

import os
import os.path as op
from itertools import product
from pprint import pprint
import pickle
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from pyplr import stlabhelp
from pysilsub.problems import SilentSubstitutionProblem as SSP
from pysilsub.CIE import get_CIE_1924_photopic_vl

LUMENS_PER_WATT = 683


# %% ~~~ PLOT STYLE ~~~

plt.style.use('seaborn')
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'Helvetica'

# %% ~~~ CONSTANTS ~~~

MINTENSITY = 0
MAXTENSITY = 4095
BACKGROUND = MAXTENSITY/2
Fs = 100  # STLAB switching time
VL = get_CIE_1924_photopic_vl()

# %% ~~~ MAKE FOLDERS FOR OUTPUT ~~~

STIM_FOLDER = './LuminanceTest/stims/'
if not op.exists(STIM_FOLDER):
    os.makedirs(STIM_FOLDER)

STLAB_1_STIM_FOLDER = op.join(STIM_FOLDER, 'STLAB_1/')
if not op.exists(STLAB_1_STIM_FOLDER):
    os.mkdir(STLAB_1_STIM_FOLDER)
    
STLAB_2_STIM_FOLDER = op.join(STIM_FOLDER, 'STLAB_2/')
if not op.exists(STLAB_2_STIM_FOLDER):
    os.mkdir(STLAB_2_STIM_FOLDER)
    
CALIBRATION_FOLDER = './calibration/'
if not op.exists(CALIBRATION_FOLDER):
    os.mkdir(CALIBRATION_FOLDER)

GAMMA_FOLDER = op.join(CALIBRATION_FOLDER, 'gamma/')
if not op.exists(GAMMA_FOLDER):
    os.mkdir(op.join(GAMMA_FOLDER))

#%% ~~~ CALIBRATION ~~~ #

# Load predictive models for each device.
S1 = SSP.from_json('./calibration/STLAB_1_York.json')
S2 = SSP.from_json('./calibration/STLAB_2_York.json')

# We know that STLAB_2 has a higher output than STLAB_1. Here we obtain
# a calibration ratio for each LED that *may* be used to perform a simple
# correction later.
S1_S2_calibration_ratio = (S1.calibration.groupby(level=0).sum().sum(axis=1)
                           / S2.calibration.groupby(level=0).sum().sum(axis=1))
print('> S1/S2 LED calibration ratio')
print(S1_S2_calibration_ratio)

# To scale Y-axis for calibration plots
max_counts = max(S1.calibration.max().max(), S2.calibration.max().max())

# Plot the calibration spds, do the gamma corrections, save output, etc.
for device in [S1, S2]:
    # Plot spds
    fig, ax = plt.subplots(figsize=(12, 4))
    device.plot_calibration_spds(ax=ax)
    ax.set_ylim(0, max_counts*1.05)
    fig.savefig(
        f'./{CALIBRATION_FOLDER}/{device.config["json_name"]}_calibration_spds.svg')

    # Keep a log of which device / calibration was used to prepare the stims
    # and at what time
    with open(f'./{CALIBRATION_FOLDER}/{device.config["json_name"]}_device_log.txt', 'w') as fh:
        pprint(S1.config, stream=fh)
        print(f'\n>> Time created: {datetime.now()}', file=fh)

    # Perform gamma correction
    device.do_gamma(fit='polynomial')
    device.gamma[device.gamma < MINTENSITY] = MINTENSITY
    device.gamma[device.gamma > MAXTENSITY] = MAXTENSITY
    device.gamma.to_csv(
        op.join(GAMMA_FOLDER, f'{device.config["json_name"]}_gamma_table.csv'))
    device.plot_gamma(save_plots_to=GAMMA_FOLDER, show_corrected=True)

# Match backgrounds and pickle
S2.background = pd.Series([.5] * S1.nprimaries)
S1.background = S2.background / S1_S2_calibration_ratio

# Pickle backgrounds so they can be loaded at start of experiment script
with open(op.join(STIM_FOLDER, 'STLAB_1_background.pickle'), 'wb') as fh:
    pickle.dump(S1.w2s(S1.background), fh)
with open(op.join(STIM_FOLDER, 'STLAB_2_background.pickle'), 'wb') as fh:
    pickle.dump(S2.w2s(S2.background), fh)

# Save plot of the background spectra
fig, ax = plt.subplots()
s1_bg = S1.predict_multiprimary_spd(S1.background)
s2_bg = S2.predict_multiprimary_spd(S2.background)
ax.plot(s1_bg, label='STLAB_1 background')
ax.plot(s2_bg, label='STLAB_2 background')
ax.legend()
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Counts/s/nm')
ax.set_title('Background spectra')
fig.savefig(op.join(STIM_FOLDER, 'Background_spectra.svg'))


#%% ~~~ MAKE THE STIMS ~~~ #

# Plot for stimuli
fig, axs = plt.subplots(2, 3, figsize=(12, 6.5))

# Stimulus parameters
frequencies = [.4, .5]
contrasts = [.06, .12, .24, .48, .96]
seconds = 12

contrast_profiles = []
for f, c in product(frequencies, contrasts):
    print(f,c)
    # A complete cycle of target contrasts
    x = stlabhelp.sinusoid_modulation(f, 1/f, Fs) * c
    t = np.linspace(0, 1/f, len(x))

    # Plot ideal stimulus profile
    if f == .5:
        axs[0, 0].plot(t, x, c='k', alpha=c, lw=2+c, label=c)
    elif f == .4:
        axs[1, 0].plot(t, x, c='k', alpha=c, lw=2+c, label=c)

    # Find the solutions
    s1_solutions = []
    s2_solutions = []
    for problem, solutions in zip([S1, S2], [s1_solutions, s2_solutions]):
        bg_spds = problem.predict_multiprimary_spd(
            problem.background, nosum=True)

        # Primary to luminance
        A = VL.T.dot(bg_spds)

        # Inverse of luminance
        A1 = 1/A

        for target in x:
            target = A * target
            # Apply the calibration ratio correction to STLAB_2
            if problem.name == 'STLAB_2 (binocular, right eye)':
                modulation = (A1 * target) * S1_S2_calibration_ratio
            else:
                modulation = (A1 * target)

            solutions.append((modulation / 2 + problem.background).squeeze())
    
    s1_solutions = [s.clip(0.0, 1.0) for s in s1_solutions]
    s2_solutions = [s.clip(0.0, 1.0) for s in s2_solutions]

            
    # Gamma correct
    s1_solutions = [S1.gamma_correct((s*(4095)).astype('int'))
                    for s in s1_solutions]
    s2_solutions = [S2.gamma_correct((s*(4095)).astype('int'))
                    for s in s2_solutions]
    
    s1ao = pd.concat([S1.get_photoreceptor_contrasts(s) for s in s1_solutions], axis=1).T
    s2ao = pd.concat([S2.get_photoreceptor_contrasts(s) for s in s2_solutions], axis=1).T
    s1ao['Device'] = 'STLAB_1'
    s2ao['Device'] = 'STLAB_2'
    s1ao['f'] = f
    s1ao['c'] = c
    s2ao['f'] = f
    s2ao['c'] = c
    ao = pd.concat([s1ao, s2ao]).reset_index().rename(columns={'index':'spectrum'})
    # Uncomment below and comment above to not use gamma
    #s1_solutions = [(s*4095).astype('int') for s in s1_solutions]
    #s2_solutions = [(s*4095).astype('int') for s in s2_solutions]

    # Put the settings in lists
    spectra1 = []
    spectra2 = []
    for s1, s2 in zip(s1_solutions, s2_solutions):
        spectra1.append(S1.predict_multiprimary_spd(s1))
        spectra2.append(S2.predict_multiprimary_spd(s2))

    # Plot forward projection of photopic luminance for all of the spectra on
    # each device
    spectra1 = pd.concat(spectra1, axis=1).T.reset_index(drop=True)
    spectra2 = pd.concat(spectra2, axis=1).T.reset_index(drop=True)
    spectra1_vl = spectra1.dot(VL) * (LUMENS_PER_WATT/100)
    spectra2_vl = spectra2.dot(VL) * (LUMENS_PER_WATT/100)
    

    if f == .5:
        axs[0, 1].plot(t, spectra1_vl.values, label=c, lw=1.5+c, alpha=c, color='k')
        axs[0, 2].plot(t, spectra2_vl.values, label=c, lw=1.5+c, alpha=c, color='k')
    elif f == .4:
        axs[1, 1].plot(t, spectra1_vl.values, label=c, lw=1.5+c, alpha=c, color='k')
        axs[1, 2].plot(t, spectra2_vl.values, label=c, lw=1.5+c, alpha=c, color='k')

    spectra1_vl['Device'] = 'STLAB_1'
    spectra2_vl['Device'] = 'STLAB_2'

    vl = pd.concat([spectra1_vl, spectra2_vl])
    vl['f'] = f
    vl['c'] = c
    vl = vl.reset_index().rename(columns={'index':'spectrum'})
    full = pd.concat([vl, ao[['sc','mc','lc','rh','mel']]], axis=1)
    contrast_profiles.append(full)

    # ~~~ MAKE STLAB VIDEO FILES ~~~ #

    time_points = np.linspace(0, 1000*seconds, seconds*Fs).astype('int')

    for device, settings in zip([S1, S2], [s1_solutions, s2_solutions]):

        settings = settings * 30  # More than we need
        settings = settings[0:len(time_points)]  # Trim it down to size
        led_cycles = np.array(settings)

        led_cycles = np.insert(led_cycles, 0, time_points, axis=1)

        # If frequency is not a round number, the device output will not
        # return to the background after 12 seconds. So we add final spectrum
        # to make sure the luminaire returns to background levels after a
        # modulation.
        final_spec = led_cycles[0].copy()
        final_spec[0] = 12050  # 50 ms afterwards, won't be skipped
        led_cycles = np.vstack([led_cycles, final_spec])

        led_cycles = pd.DataFrame(led_cycles)
        led_cycles.columns = ['time' if c == 0
                              else 'LED-' + str(c-1)
                              for c in led_cycles.columns]
        metadata = {
            'title': f'{f} Hz luminance modulation',
            'seconds': seconds,
            'contrast': c,
            'frequency': f,
            'device': device.name
        }
        fname = op.join(
            STIM_FOLDER, f'{device.name.split()[0]}/c{c}_f{f}.json')
        stlabhelp.make_video_file(
            led_cycles, repeats=1, fname=fname, **metadata)

# Tweak the stim figure and save
for ax in [axs[0, 1], axs[0, 2], axs[1, 1], axs[1, 2]]:
    ax.set_ylabel('Illuminance (lux)')
    ax.set_ylim((0, 160))

axs[0, 0].set_ylabel('Target luminance contrast\n$f$ = .5')
axs[1, 0].set_ylabel('Target luminance contrast\n$f$ = .4')
axs[0, 0].set_title('Ideal stimulus profile')
axs[0, 1].set_title('Forward prediction - Source 1')
axs[0, 2].set_title('Forward prediction - Source 2')
axs[0, 2].legend(title='Contrast')

for ax in axs[1]:
    ax.set_xlabel('Time (s)')

plt.tight_layout()

fig.savefig(op.join(STIM_FOLDER, 'STLAB_luminance_contrast_stims.svg'))
contrast_profiles = pd.concat(contrast_profiles, axis=0)
contrast_profiles.to_csv(op.join(STIM_FOLDER, 'contrast_profiles.csv'))
