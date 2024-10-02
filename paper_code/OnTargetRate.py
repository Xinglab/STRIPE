#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.02

Compute median on-target rate for both CDG and PMD gene panels for TEQUILA-seq samples in manuscript.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import numpy as np
import pandas as pd

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = "/mnt/isilon/lin_lab_share/STRIPE"

# =====================================================================================================================
#                                                     CDG PANEL
# =====================================================================================================================

sampleDF, rates = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0), []

for sampid in sampleDF[sampleDF['Provider'] != 'Lan Lin']['ID'].values:
    rates.append(float(pd.read_csv(workdir + '/CDG/' + sampid + '/RNA/stripe/quality_control/' + sampid + 
        '_TEQUILA.mapping_stats.txt', sep = '\t', header = None).tail(1)[0].item().replace('%', '')))

print('Median on-target rate (CDG gene panel): ' + str(np.median(rates)) + '%')

# =====================================================================================================================
#                                                     PMD PANEL
# =====================================================================================================================

sampleDF, rates = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0), []

for sampid in sampleDF[sampleDF['Provider'] != 'Rebecca Ganetzky']['ID'].values:
    rates.append(float(pd.read_csv(workdir + '/PMD/' + sampid + '/RNA/stripe/quality_control/' + sampid + 
        '_TEQUILA.mapping_stats.txt', sep = '\t', header = None).tail(1)[0].item().replace('%', '')))

print('Median on-target rate (PMD gene panel): ' + str(np.median(rates)) + '%')
