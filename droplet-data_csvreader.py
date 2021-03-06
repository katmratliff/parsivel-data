# -*- coding: utf-8 -*-
"""
PURPOSE: - import raw OTT Parsivel disdrometer data
         - exports .csv file with precip statistics
NOTES:   - USER MUST UPDATE FILENAME & RESAVE before running
         - Raw data file (*.MIS) must be in same folder as this script
         - must rename .MIS files to remove periods before extension
         - user can specify number of rows of data to be deleted from stats
         - user can specify volume collected water to calculate % error
         - user can specify tolerance range for excluding outliers
author: Katherine Ratliff (ratliff.katherine@epa.gov)
"""

import pandas as pd
import numpy as np
import csv
import pdb

FILENAME = 'Protexusbp2.MIS'

# amount of water (in mL) collected in bin during test
# and bin dimensions (used to calculate % error)
# and collection time for water (min)
vol_collected = 0
bin_dim = 9.75*12
coll_time = 10

# delete first row(s) of data from MIS file from calcs?
# number specifies how many rows to remove
del_top_rows = 2

# delete last row(s) of data from MIS file for calcs?
# number specifies how many rows to remove
del_bot_rows = 0

# what is your outlier tolerance? i.e., how many standard deviations away
# from the mean do we keep data (to avoid erroneous parsivel data)
# if zscore_tol = 0, no data excluded
zscore_tol = 5

mm_to_in = 0.0393701 # converts mm to inches

particle_sizes = np.loadtxt('particle-sizes.txt')

df = pd.DataFrame()

# import and split dataframe
#df = pd.read_csv(FILENAME, header=None, index_col=False, engine='python', sep=',',
#                 skiprows = lambda x: x in [0, del_top_rows+1])

with open (FILENAME, newline='') as csvfile:
    misreader = csv.reader(csvfile, delimiter=',') #, quotechar='|',
                           #quoting=csv.QUOTE_MINIMAL)
    #pdb.set_trace()
    for row in misreader:
        if row == []:
            pass
        else:
            ds = pd.Series(row)
            df = df.append(ds, ignore_index=True)
            #df = pd.concat([df, ds], ignore_index=True)
        #print(row)

# delete number of top and/or bottom rows if specified
if del_top_rows or del_bot_rows:
    df = df.iloc[del_top_rows:(len(df)-del_bot_rows)]

# delete remaining rows where spectral data doesn't exist
df = df.drop(df[df.iloc[:, 8] != '<SPECTRUM>'].index)

# split out droplet data and add column headers
drop_data = df.iloc[:, :8]
drop_data.columns= ["date", "time", "intensity", "total_precip",
                    "reflectivity", "visibility", "num_particles", "KE"]

# split out spectral data and fill with zeros as needed
spec_raw = df.iloc[:, 8:1032]
spec_raw = spec_raw.replace('<SPECTRUM>', 0)
spec_raw = spec_raw.fillna(0)

# initialize spectral data array
spectrum = np.zeros([len(spec_raw), 32, 32])

# reshape raw spectrum data
for i in range(0, len(spec_raw)):
    spectrum[i] = spec_raw.iloc[i, :].values.reshape(32, 32)

# sum number of droplets per size class for each measurement
num_drops_per_diam = np.sum(spectrum, axis=1)

# multiply number of dropelets per size class by its size
weighted_size = num_drops_per_diam * particle_sizes

# calculate average diameter for each period
avg_diam_overtime = np.divide(np.sum(weighted_size, axis=1),
                              np.sum(num_drops_per_diam, axis=1))

# calculate average droplet size weighted by diameter over dataset
# (standard deviation is variation in average diameter over time)
weighted_avg_diam = np.divide(np.sum(avg_diam_overtime), len(spec_raw))
avg_diam_std = np.std(avg_diam_overtime)

# remove rows where data falls outside tolerated z-score
# for now, only looking at number of particles to do this
# (can add more/different criteria later)
if zscore_tol:
    drop_data = drop_data[~(np.abs(drop_data.num_particles -
                                  drop_data.num_particles.mean())
                          >(zscore_tol * drop_data.num_particles.std()))]

# stats from parsivel
intensity_avg_mmhr = np.mean(drop_data.intensity)
intensity_std_mmhr = np.std(drop_data.intensity, ddof=1)

intensity_avg_inhr = np.mean(drop_data.intensity * mm_to_in) 
intensity_std_inhr = np.std(drop_data.intensity * mm_to_in, ddof=1)

num_particles_avg = np.mean(drop_data.num_particles)
num_particles_std = np.std(drop_data.num_particles, ddof=1)

kinetic_energy_avg = np.mean(drop_data.KE)
kinetic_energy_std = np.std(drop_data.KE, ddof=1)

# collected water stats, conversion done as in excel sheet
if vol_collected:
    collected_intensity = (vol_collected * 0.0610237 * 60 / (coll_time*bin_dim))
    perc_error = np.abs(1 - collected_intensity / intensity_avg_inhr)

# write ouptut to rainfall_stats.csv
with open('rainfall_stats.csv', 'w', newline='') as csvfile:
    precip = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    precip.writerow(['Parsivel Statistics'])
    precip.writerow(['Variable', 'Average', 'Standard Deviation'])
    precip.writerow(['Rain Intensity [mm/hr]', intensity_avg_mmhr,
                     intensity_std_mmhr])
    precip.writerow(['Rain Intensity [in/hr]', intensity_avg_inhr,
                     intensity_std_inhr])
    precip.writerow(['Number of Particles', num_particles_avg,
                     num_particles_std])
    precip.writerow(['Weighted Average Diameter [mm]', weighted_avg_diam,
                     avg_diam_std])
    precip.writerow(['Kinetic Energy', kinetic_energy_avg, kinetic_energy_std])
    if vol_collected:
        precip.writerow([''])
        precip.writerow(['Collected Water Statistics'])
        precip.writerow(['Volume (mL)', vol_collected])
        precip.writerow(['Rain Intensity (in/hr)', collected_intensity])
        precip.writerow(['Percent Error', "{:.2%}".format(perc_error)])
