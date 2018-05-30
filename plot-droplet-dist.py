# -*- coding: utf-8 -*-
"""
PURPOSE: - import raw OTT Parsivel disdrometer data
         - plots particle size, velocity, & volumetric distributions
NOTEs:   - USER MUST UPDATE FILENAME/NOZZLE & RESAVE before running
         - Raw data file (*.MIS), 'particle-sizes.txt' & 'velocity-bins.txt'
           must be in same folder as this script
         - plots are placed in a 'figures' subfolder
         - plot axes ranges can be specified below
author: Katherine Ratliff (ratliff.katherine@epa.gov)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

FILENAME = 'FL5-VS, 50 mesh, 50 psi 1 hr test.MIS'
nozzle = 'TG-1'

# delete first row of data from MIS file from calcs? 1 = yes, 0 = no
del_row = 0

# plot maximum std deviation on figures? (for size/velocity)
plot_errorbars = 0

particle_sizes = np.loadtxt('particle-sizes.txt')
velocity_bins = np.loadtxt('velocity-bins.txt')

if not os.path.exists("figures"):
    os.mkdir("figures")

# import and split dataframes
df = pd.read_csv(FILENAME, header=None)
if del_row:
    df_meta = df.iloc[1:, :8]
    df_spec_raw = df.iloc[1:, 8:1032]
else:
    df_meta = df.iloc[:, :8]
    df_spec_raw = df.iloc[:, 8:1032]
df_spec_raw = df_spec_raw.replace('<SPECTRUM>', 0)
df_spec_raw = df_spec_raw.fillna(0)

# initialize spectral data arrays
spectrum = np.zeros([len(df_spec_raw), 32, 32])
spec_avg = np.zeros([32, 32])

for i in range(0, len(df_spec_raw)):
    
    # reshape raw spectrum data
    spectrum[i] = df_spec_raw.iloc[i, :].values.reshape(32, 32)

# average and std. dev. of spectrum data over time
spec_avg = np.average(spectrum, axis=0)
spec_std = np.std(spectrum, axis=0)

# sum spectral data over time to use in calculating intensity
spec_sum = np.sum(spectrum, axis=0)

# sum spectral data over velocity axis
size_sum = np.sum(spec_avg, axis=0)
size_std_max = np.amax(spec_std, axis=0)

# total number of droplets per diameter
tot_drops = np.sum(spec_sum, axis=0)

# convert droplet count to volume
vol_per_diam = tot_drops * ((4/3) * np.pi * np.power((particle_sizes/2), 3))
tot_volume = np.sum(vol_per_diam)
vol_norm = vol_per_diam / tot_volume

# calculate volume mean diameter (VMD)
cum_vol = np.cumsum(vol_per_diam)
VMD_loc = len(cum_vol[cum_vol < (tot_volume/2)])
VMD = particle_sizes[VMD_loc]

# sum spectral data over size axis
velocity_sum = np.sum(spec_avg, axis=1)
velocity_std_max = np.amax(spec_std, axis=1)

# remove zeros from data so they don't plot
zeros_index_size = np.where(size_sum == 0)
size_sum = np.delete(size_sum, zeros_index_size)
particle_sizes = np.delete(particle_sizes, zeros_index_size)
size_std_max = np.delete(size_std_max, zeros_index_size)
vol_norm = np.delete(vol_norm, zeros_index_size)

zeros_index_velocity = np.where(velocity_sum == 0)
velocity_sum = np.delete(velocity_sum, zeros_index_velocity)
velocity_bins = np.delete(velocity_bins, zeros_index_velocity)
velocity_std_max = np.delete(velocity_std_max, zeros_index_velocity)

if plot_errorbars: 
    # plot size distribution figure
    f = plt.figure()
    plt.errorbar(particle_sizes, size_sum, size_std_max, fmt='o')
    plt.ylim((-50, 1300))
    plt.xscale('log')
    plt.xlabel('log of droplet diameter (mm)')
    plt.ylabel('number of droplets')
    plt.title(FILENAME[:-4])
    plt.savefig('figures/' + nozzle + '_size_dist.png')
    plt.close(f)
    
    # plot velocity distribution figure
    f = plt.figure()
    plt.errorbar(velocity_bins, velocity_sum, velocity_std_max,
                 fmt='o', color='r')
    plt.ylim((-20, 800))
    plt.xscale('log')
    plt.xlabel('log of particle velocity (m/s)')
    plt.ylabel('number of droplets')
    plt.title(FILENAME[:-4])
    plt.savefig('figures/' + nozzle + '_velocity_dist.png')
    plt.close(f)
    
else:
    # plot size distribution figure
    f = plt.figure()
    plt.scatter(particle_sizes, size_sum)
    plt.ylim((-50, 1300))
    plt.xscale('log')
    plt.xlabel('log of droplet diameter (mm)')
    plt.ylabel('number of droplets')
    plt.title(FILENAME[:-4])
    plt.savefig('figures/' + nozzle + '_size_dist.png')
    plt.close(f)
    
    # plot velocity distribution figure
    f = plt.figure()
    plt.scatter(velocity_bins, velocity_sum, color='r')
    plt.ylim((-20, 800))
    plt.xscale('log')
    plt.xlabel('log of particle velocity (m/s)')
    plt.ylabel('number of droplets')
    plt.title(FILENAME[:-4])
    plt.savefig('figures/' + nozzle + '_velocity_dist.png')
    plt.close(f)
    
# plot normalized volume figure
f = plt.figure()
plt.scatter(particle_sizes, vol_norm)
plt.text(0.3, 0.5, ('VMD = ' + str(np.max(VMD)) + ' mm'))
plt.axis([0.25, 25, -0.05, 0.6])
plt.xscale('log')
plt.xlabel('log of droplet diameter (mm)')
plt.ylabel('fraction of total rainfall volume')
plt.title(FILENAME[:-4])
plt.savefig('figures/' + nozzle + '_volume.png')
plt.close(f)

