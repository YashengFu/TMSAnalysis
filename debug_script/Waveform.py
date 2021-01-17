import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cycler

import time
import sys
sys.path.append('/dybfs2/nEXO/fuys/stanford_teststand')
import os

from TMSAnalysis.StruckAnalysisConfiguration import StruckAnalysisConfiguration
from TMSAnalysis.WaveformAnalysis import Waveform

plt.rcParams['axes.prop_cycle'] = cycler.cycler(color='bgrmyk')
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
plt.rcParams['figure.figsize'] = [10, 8]
plt.rcParams['font.size'] = 12

#run_parameters_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Run_Parameters_Xe_Run30_SimCompatible.csv'
run_parameters_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Run_Parameters_Run30_20200915_Night_AfterFirstRnInjection.csv'
calibrations_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Calibrations_Xe_Run30.csv'
#channel_map_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Channel_Map_Run30.csv'
channel_map_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/31th/Channel_Maps_Run31_DS01.csv'
#channel_map_file = '/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Channel_Map_Run30_20200916_RnPoAlphaEffTest.csv'


#analysis_config object loads all these paramenters
analysis_config = StruckAnalysisConfiguration.StruckAnalysisConfiguration()
analysis_config.GetRunParametersFromFile( run_parameters_file )
analysis_config.GetChannelMapFromFile( channel_map_file )

#for parameter, value in analysis_config.run_parameters.items():
#    print( '{:<25}\t{:>10.6}'.format( parameter+':', float(value) ) )
    
sampling_period_us = analysis_config.run_parameters['Sampling Period [ns]'] / 1.e3 # Convert from ns to us
trigger_time_samples = analysis_config.run_parameters['Pretrigger Length [samples]']
waveform_length_samples = analysis_config.run_parameters['Waveform Length [samples]']

#path_to_reduced_data = '/dybfs2/nEXO/fuys/stanford_teststand/data/30th/20200912_MorningNoise_PreRecirculation/analysis_500ns/'
path_to_reduced_data = '/dybfs2/nEXO/fuys/stanford_teststand/data/31th/20201102_DS23_NoiseData/analysis_500ns/'

#fname = path_to_reduced_data + 'tier1_SIS3316Raw_20200912180038_SiPMs_longTPC_sbias33p0_scope_trig13_35mV_cath_6022V_1-ngm_reduced.h5'
fname = path_to_reduced_data + 'tier1_SIS3316Raw_20201102165346Run31_DS23_SiPMs_longTPC_sbias32p5_internalTrigger150ADC_3foldCoin_cath_6000V_Noise_1-ngm_reduced.h5'

data_df = pd.read_hdf(fname)
#data_df.head()

#for colname in data_df.columns:
#    print(colname)

#charge_energy = data_df['TotalTileEnergy']
#time_of_max_channel = data_df['TimeOfMaxChannel']

#data can also be selected for example a specific energy range
#mask = (data_df['TotalTileEnergy']>100) & (data_df['TotalTileEnergy']<5000)& (data_df['TotalSiPMEnergy']>100) & (data_df['TotalSiPMEnergy']<40000)
#make a 2d histogram with the selected data (for example charge energy vs light energy)

#xbins = np.linspace(0.,6500.,200)
#ybins = np.linspace(0.,5000.,100)

#plt.hist2d(time_of_max_channel[mask], charge_energy[mask],bins=(xbins,ybins))
#plt.plot(time_of_max_channel[mask], charge_energy[mask],'o',color=(0.,0.,1.,0.1),markersize=1.)
#plt.xlabel('Time of max channel (samples)')
#plt.ylabel('Charge Energy (ADC units)')

#drift_time = (time_of_max_channel - analysis_config.run_parameters['Pretrigger Length [samples]'])*analysis_config.run_parameters['Sampling Period [ns]'] / 1000.

#plt.hist2d(time_of_max_channel[mask], charge_energy[mask],bins=(xbins,ybins))
#plt.plot(drift_time[mask], charge_energy[mask],'o',color=(0.,0.,1.,0.1),markersize=1.)
#
#plt.xlabel('Drift time (microseconds)')
#plt.ylabel('Charge Energy (ADC units)')

# Next, we can select events where the charge is only on one or two channels. This should be the case for most alpha events, since the ionization is so localized. 

mask = (data_df['TotalTileEnergy']>0) & (data_df['TotalTileEnergy']<500000)&(data_df['TotalSiPMEnergy']<100000) & (data_df['TotalSiPMEnergy']<4000000)&(data_df['NumTileChannelsHit'] < 3000)

#plt.hist2d(time_of_max_channel[mask], charge_energy[mask],bins=(xbins,ybins))
#plt.plot(drift_time[mask], charge_energy[mask],'o',color=(0.,0.,1.,0.1),markersize=1.)
#
#plt.xlabel('Drift time (microseconds)')
#plt.ylabel('Charge Energy (ADC units)')

#the indexes of the events passing the cuts can be seen in this way 
cut_index = data_df[mask].index
#print(data_df[mask])
#print('----cut_index[0]------')
print(cut_index)
#print('-----mask----')
print(mask)
#a waveform can be accessed in the following way, in this case the first event passing the cut
event = Waveform.Event(fname,\
    #'/dybfs2/nEXO/fuys/stanford_teststand/data/30th/20200912_MorningNoise_PreRecirculation/raw_data',\
    '/dybfs2/nEXO/fuys/stanford_teststand/data/31th/20201102_DS23_NoiseData/raw_data',\
    cut_index[30],\
    run_parameters_file,calibrations_file,channel_map_file)
#the second parameter in the Event class is the location of the tier1 file, which is where the low level information
#are stored, such as the waveform. Make sure the dataset (in this case 20200213_AfterRnInjection matches with the one
#opened in the reduced_added.h5 file), the third parameter is the number of the event you want to see the waveform of

smoothing_windows_us = 0.08
#event.smooth(smoothing_windows_us)
#event.FFT()
plot = event.plot_event(risetime=True)
#plot.show()
plot.savefig('./waveform.png')
