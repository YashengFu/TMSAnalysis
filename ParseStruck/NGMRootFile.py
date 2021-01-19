############################################################################
# This file defines a class that reads the ROOT files produced by the 
# NGM Daq and puts the data into Pandas dataframes 
# for analysis purposes.
#
#    - Brian L.
############################################################################

import pandas as pd
import numpy as np
import uproot3 as up
import time
import os
import sys
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

class NGMRootFile:

        ####################################################################
	def __init__( self, input_filename=None, output_directory=None, config=None, start_stop = [None, None]):
		print('NGMFile object constructed.')

		self.start_stop = start_stop
		package_directory = os.path.dirname(os.path.abspath(__file__))
		if output_directory is None:
			self.output_directory = './'
		else:
			self.output_directory = output_directory + '/'

		if input_filename is not None:
			self.LoadRootFile( input_filename )
		if config is not None:
			self.channel_map = config.channel_map
		else:
			print('WARNING: No channel map file provided. Using the default one...')
			#self.channel_map = pd.read_csv('/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Channel_Map_Run30.csv',skiprows=0)
			self.channel_map = pd.read_csv('/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/31th/Channel_Maps_Run31_DS01.csv',skiprows=0)
			#self.channel_map = pd.read_csv(package_directory + '/channel_map_8ns_sampling.txt',skiprows=9)
		self.h5_file = None

        ####################################################################
	def LoadRootFile( self, filename ):
		self.infile = up.open(filename)
		self.filename = filename
		print('Input file: {}'.format(self.filename))
		try:
			self.intree = self.infile['HitTree']
		except ValueError as e:
			print('Some problem getting the HitTree out of the file.')
			print('{}'.format(e))
			return
		print('Got HitTree.')
		

        ####################################################################
	def GroupEventsAndWriteToHDF5( self, nevents = -1, save = True, start_stop = None ):
		
		try:
			self.infile
		except NameError:
			self.LoadFile()

		if start_stop is not None:
			self.start_stop = start_stop	

		start_time = time.time()		
		file_counter = 0
		global_evt_counter = 0
		local_evt_counter = 0
		df = pd.DataFrame(columns=['Channels','Timestamp','Data','ChannelTypes','ChannelPositions','SiPMEnergy'])
		start_time = time.time()
		print('{} entries per event.'.format(len(self.channel_map)))

		for data in self.intree.iterate(['_waveform','_rawclock','_slot','_channel'],\
						namedecode='utf-8',\
						entrysteps=len(self.channel_map),\
						entrystart=self.start_stop[0],\
						entrystop=self.start_stop[1]):
			if nevents > 0:
				if global_evt_counter > nevents:
					break

			data_series = pd.Series(data)
			channel_mask, channel_types, channel_positions = self.GenerateChannelMask( data['_slot'],data['_channel'])
                         
                        # Remove 'Off' channels from the data stream
			for column in data_series.items():
				data_series[ column[0] ] = np.array(data_series[column[0]][channel_mask])
			output_series = pd.Series()
			output_series['Channels'] = data_series['_slot']*16+data_series['_channel']
			output_series['Timestamp'] = data_series['_rawclock']
			output_series['Data'] = data_series['_waveform']#by yasheng
			channel_mask, channel_types, channel_positions = self.GenerateChannelMask( data_series['_slot'],data_series['_channel'])
			output_series['ChannelTypes'] = channel_types
			#print('-----------------------')
			#print(type(data_series['_waveform'][0]))
			#print(type(output_series['Channels']))
			#print(len(output_series['Channels']))
			output_series['ChannelPositions'] = channel_positions
			###################################################################
		    #for SiPM gain calibration Yasheng Fu 12/23/2020
			'''
			sampling_period_ns = 1./(62.5/1.e3)#hard code need to change
			datas=[]
			for i in np.arange(len(output_series['Channels'])):
				peak_energy=[]
				if 'SiPM' in output_series['ChannelTypes'][i]:# and output_series['Channels'] == 0:
					baseline_value = np.mean(data_series['_waveform'][i])	
					#print(baseline_value)
					#Transform to frequency domain
					wfm_fft =np.fft.rfft(data_series['_waveform'][i]-baseline_value)
					
					#FFT and high frequence cut(<0.5e7)
					wfm_fft_pass = np.zeros_like(wfm_fft) 
					freqs = np.fft.rfftfreq(len(data_series['_waveform'][i]),d=1.e-9*sampling_period_ns)
					fft_freq_pass = np.logical_and(freqs < 0.5e7, freqs >=0)
					wfm_fft_pass[fft_freq_pass] = wfm_fft[fft_freq_pass]	

					#inverse FFT
					wfm_ifft_pass = np.fft.irfft(wfm_fft_pass)

					#Peak finding algorithm
					peaks, _= find_peaks(wfm_ifft_pass,height = 150)
					#print("---------------peaks-------------")
					#print(type(peaks))
					if len(peaks):
						for i in np.arange(len(peaks)):
							t_max_point = peaks[i]
							window_start = t_max_point - int(200/sampling_period_ns)
							window_end = t_max_point + int(200/sampling_period_ns)
							data_array = wfm_ifft_pass[window_start:window_end]
							cumul_pulse_energy = np.cumsum(data_array)
							area_window_length = int(5./sampling_period_ns)# average over 5ns
							sipm_energy = np.mean(cumul_pulse_energy[-area_window_length:])
							#peak_energy.append(sipm_energy)
							peak_energy.append(wfm_ifft_pass[peaks[i]])#pluse test by yasheng
				datas.append(peak_energy)
						#print("Not empty")ls
			'''
					#print(peaks) 
			'''
				data_series['_waveform'][i] = gaussian_filter( data_series['_waveform'][i].astype(float), \
				80./sampling_period_ns )
				baseline = np.mean(data_series['_waveform'][i])
				max_point = np.max(data_series['_waveform'][i])
				t_max_point = np.argmax(data_series['_waveform'][i])
				
				window_start = t_max_point - int(400/sampling_period_ns)
				window_end = t_max_point + int(450/sampling_period_ns)

				data_array = data_series['_waveform'][i][window_start:window_end]-baseline
				cumul_pulse_energy = np.cumsum( data_array)
				area_window_length = int(50./sampling_period_ns)# average over 50ns
				sipm_energy = np.mean(cumul_pulse_energy[-area_window_length:])
				#print(max_point)
				#print(max_point-baseline)
				#print(sipm_energy) 
				datas.append(sipm_energy)
			output_series['SiPMEnergy'] = datas			#print(type(output_series['SiPMEnergy']))
			'''
            ###################################################################
			df = df.append(output_series,ignore_index=True)	


			global_evt_counter += 1
			local_evt_counter += 1
			#'''by yasheng
			if local_evt_counter > 100 and save:
				output_filename = '{}{}_{:0>3}.h5'.format( self.output_directory,\
									self.GetFileTitle(str(self.infile.name)),\
									file_counter)
				df.to_hdf(output_filename,key='raw')
				df.drop(df.index,inplace=True)
				local_evt_counter = 0
				file_counter += 1
				#df = pd.DataFrame(columns=['Channels','Timestamp','Data','ChannelTypes','ChannelPositions','SiPMEnergy'])
				print('Written to {} at {:4.4} seconds'.format(output_filename,time.time()-start_time))	
			#'''by yasheng

		if save:
			output_filename = '{}{}_{:0>3}.h5'.format( self.output_directory,\
								self.GetFileTitle(str(self.infile.name)),\
								file_counter )
			df.to_hdf(output_filename,key='raw')
			end_time = time.time()
			print('{} events written in {:4.4} seconds.'.format(global_evt_counter,end_time-start_time))
		else:
			return df
	
	####################################################################
	def GenerateChannelMask( self, slot_column, channel_column ):

		channel_mask = np.array(np.ones(len(slot_column),dtype=bool))
		channel_types = ['' for i in range(len(slot_column))]
		channel_positions = np.zeros(len(slot_column),dtype=int)

		for index,row in self.channel_map.iterrows():
			
			slot_mask = np.where(slot_column==row['Board'])
			chan_mask = np.where(channel_column==row['InputChannel'])
			intersection = np.intersect1d(slot_mask,chan_mask)
			if len(intersection) == 1:
				this_index = intersection[0]
			else:
				 continue
			channel_types[this_index] = row['ChannelType']
			channel_positions[this_index] = row['ChannelPosX'] if row['ChannelPosX'] != 0 else row['ChannelPosY']
			if row['ChannelType']=='Off':
				channel_mask[this_index] = False
		return channel_mask, channel_types, channel_positions
	
        ####################################################################
	def GetFileTitle( self, filepath ):
		filename = filepath.split('/')[-1]
		filetitle = filename.split('.')[0]
		return filetitle
		
	####################################################################
	def GetTotalEntries( self ):
		return self.intree.numentries
