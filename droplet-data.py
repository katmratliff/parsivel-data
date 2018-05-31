# -*- coding: utf-8 -*-
"""
PURPOSE: - import raw OTT Parsivel disdrometer data
         - exports .csv file with precip statistics
NOTES:   - USER MUST UPDATE FILENAME & RESAVE before running
         - Raw data file (*.MIS) must be in same folder as this script
         - user can specify number of rows of data to be deleted from stats
         - user can specify volume collected water to calculate % error
         - user can specify tolerance range for excluding outliers
author: Katherine Ratliff (ratliff.katherine@epa.gov)
"""

import pandas as pd
import numpy as np
import csv

FILENAME = 'TEST FILE.MIS'

# amount of water (in mL) collected in bin during test
# (used to calculate % error)
vol_collected = 0

# delete first row of data from MIS file from calcs? 1 = yes, 0 = no
del_rows = 0

# what is your outlier tolerance? i.e., how many standard deviations away
# from the mean do we keep data (to avoid erroneous parsivel data)
# if zscore_tol = 0, no data excluded
zscore_tol = 5

mm_to_in = 0.0393701 # converts mm to inches

# import and split dataframe
df = pd.read_csv(FILENAME, header=None)
if del_rows:
    drop_data = df.iloc[del_rows:, :8]
else:
    drop_data = df.iloc[:, :8]

drop_data.columns= ["date", "time", "intensity", "total_precip",
                    "reflectivity", "visibility", "num_particles", "KE"]

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
    collected_intensity = (vol_collected * 0.0610237 * 60 / (20*12.25*12.25))
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
    precip.writerow(['Kinetic Energy', kinetic_energy_avg, kinetic_energy_std])
    if vol_collected:
        precip.writerow([''])
        precip.writerow(['Collected Water Statistics'])
        precip.writerow(['Volume (mL)', vol_collected])
        precip.writerow(['Rain Intensity (in/hr)', collected_intensity])
        precip.writerow(['% Error', "{:.2%}".format(perc_error)])
