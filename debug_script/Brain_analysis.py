import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import histlite as hl
plt.rcParams['figure.figsize'] = [10, 10]
plt.rcParams.update({'font.size': 18}) 
from scipy.ndimage import gaussian_filter
import cycler
import sys
sys.path.append('/dybfs2/nEXO/fuys/stanford_teststand/')
from TMSAnalysis.StruckAnalysisConfiguration import StruckAnalysisConfiguration
import os
import time
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

analysis_config = StruckAnalysisConfiguration.StruckAnalysisConfiguration()
analysis_config.GetChannelMapFromFile('/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Channel_Map_Run30.csv')
analysis_config.GetRunParametersFromFile('/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Run_Parameters_Run30_20200911_OvernightBi207_AfterFilling.csv')
analysis_config.GetCalibrationConstantsFromFile('/dybfs2/nEXO/fuys/stanford_teststand/TMSAnalysis/config/30th/Calibrations_Xe_Run11b.csv')

#reading HDF5 filelist
files = os.listdir('/dybfs2/nEXO/fuys/stanford_teststand/data/30th/20200912_MorningNoise_PreRecirculation/Cali_data/data/')
#print(len(files))
dflist = []
file_count =0 
for filename in files:
    #if file_count%10==0:
    print(file_count)
    if file_count > 50:#by yasheng
        break
    file_count +=1
    if '0912' in filename and filename.endswith('.h5'): #and file_count <50:
        dflist.append(pd.read_hdf('/dybfs2/nEXO/fuys/stanford_teststand/data/30th/20200912_MorningNoise_PreRecirculation/Cali_data/data/'+filename))
        
df05 = pd.concat(dflist,ignore_index=True)

def Gaussian(x,A,mu,sig):
    return A*np.exp(-(x-mu)**2/(2*sig**2))


#display raw data waveform
evt_num = 3
evt = df05.iloc[evt_num]

plt.rcParams['figure.figsize'] = [12,12]
plt.rcParams.update({'font.size': 18})

chlabels = []
yvalues = []
num_channels = 0

for i in range(len(evt['Channels'])):
    if 'SiPM' not in evt['ChannelTypes'][i]:
        continue
    else:
        num_channels += 1
    
    
    ch_name = analysis_config.GetChannelNameForSoftwareChannel(evt['Channels'][i])
    
    chlabels.append(ch_name)
    yvalues.append((num_channels-1)*500.)
    
    baseline = np.mean(evt['Data'][i][0:200]) 
    
    print(evt['Data'][i])
    

    plt.plot(evt['Data'][i]-baseline + yvalues[num_channels-1],\
             label='Ch {}'.format(ch_name),color=color_cycle[num_channels-1])
    
    
plt.xlim(-100.,12000.)    
plt.legend(fontsize=18,loc='upper right',handlelength=0.75)
plt.xlabel('Time (samples)')
plt.ylabel('Voltage ADC counts (500ADC offset)')
plt.title('Raw waveforms')
plt.savefig('./Raw_data.png')
plt.close()


#display waveforms after process
from scipy.signal import butter
from scipy.signal import filtfilt
import scipy.optimize as opt
#evt_num =0
evt = df05.iloc[evt_num]
print(evt_num)

plt.rcParams['figure.figsize'] = [12,12]
plt.rcParams.update({'font.size': 18})

chlabels = []
yvalues = []
num_channels = 0

for i in range(len(evt['Channels'])):
    if 'SiPM' not in evt['ChannelTypes'][i]:
        continue
    else:
        num_channels += 1   
    
    ch_name = analysis_config.GetChannelNameForSoftwareChannel(evt['Channels'][i])
    
    chlabels.append(ch_name)
    yvalues.append((num_channels-1)*500.)
    
    baseline = np.mean(evt['Data'][i][0:200])   
    raw_data = evt['Data'][i]-baseline  
    nyq = 0.5 * 1./len(raw_data)  
    
    filter_params = butter(5,0.03,btype='highpass')
    smoothed_data = filtfilt(filter_params[0],filter_params[1],raw_data)   
    smoothed_data = gaussian_filter(smoothed_data,2.0)

    plt.plot(smoothed_data + yvalues[num_channels-1],\
             label='Ch {}'.format(ch_name),color=color_cycle[num_channels-1])
       
plt.xlim(-100.,12000.)    
plt.legend(fontsize=18,loc='upper right',handlelength=0.75)
plt.xlabel('Time (samples)')
plt.ylabel('Voltage ADC counts (500ADC offset)')
plt.title('after filter and smooth')
plt.savefig('./one_event_after_process.png')
plt.close()
#process data, include butter filter and smooth 
adc_vals = dict()
print(chlabels)
for label in chlabels:
    adc_vals[label] = []
    
for index, evt in df05.iterrows():
    if index % 100 == 0:
        print('Running event {}'.format(index))
    for i in range(len(evt['Channels'])):
        if 'SiPM' not in evt['ChannelTypes'][i]:
            continue
        else:
            num_channels += 1
        row = (num_channels-1) % 3
        col = int(np.floor( (num_channels-1)/3. ))

        ch_name = analysis_config.GetChannelNameForSoftwareChannel(evt['Channels'][i])

        baseline = np.mean(evt['Data'][i][0:200])
        raw_data = evt['Data'][i]-baseline        
        
        filter_params = butter(5,0.03,btype='highpass')
        smoothed_data = filtfilt(filter_params[0],filter_params[1],raw_data)
        smoothed_data = gaussian_filter(smoothed_data,2.0)
        
        adc_vals[ch_name].extend( list(smoothed_data) )
        
for channel, data in adc_vals.items():
    adc_vals[channel] = np.array(data)       
print(adc_vals)


#fix sigma value for find peak algorithm(use smoothed data)
sigmas_smoothed = dict()

plt.rcParams['figure.figsize'] = [18,14]
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(ncols=4,nrows=3,sharex=True,sharey=True,gridspec_kw={'hspace':0.,'wspace':0.})

index = 0
for channel, array in adc_vals.items():

    row = (index) % 3
    col = int(np.floor( index/3. ))
    
    thishist = hl.hist( array, bins=np.linspace(-300.,1000.,301) )
    hl.plot1d( ax[row,col], thishist, color=color_cycle[index],label=channel )
    
    ax[row,col].legend(fontsize=14)
    ax[row,col].set_yscale('log')
    ax[row,col].set_ylim(0.5,5000000.)
    ax[row,col].set_xlim(-300.,1000.)
    
    bin_centers = (thishist.bins[0][1:] + thishist.bins[0][:-1])/2.
    bin_vals = thishist.values
    fitmask = (bin_centers>-100.)&(bin_centers<100.)

    if '2-1'  not in channel and '2-3' not in channel and '2-4' not in channel:
    #if '1-3' not in channel and '1-1' not in channel:
        p,pcov = opt.curve_fit(Gaussian,bin_centers[fitmask],bin_vals[fitmask],p0=(5000000.,0.,10.))
        sigmas_smoothed[channel] = p[2]
        print('Channel {}: 1sigma = {:4.4}'.format(channel,p[2]))
        textstring = r'$\sigma$ = {:3.3} ADC'.format(p[2])
        ax[row,col].text(70.,200.,textstring,fontsize=12)

        xfit = np.linspace(-300.,300.,200)
        yfit = Gaussian(xfit,p[0],p[1],p[2])
        ax[row,col].plot(xfit,yfit,'--k',linewidth=1)
    
    index += 1

plt.savefig('./histogram_of_sigma_display_5000events.png',dpi=200,bbox_inches='tight')
plt.close()

# define function to find singal in noise waveforms
def PulseFinder( raw_data, chname, threshold_sig=4. ):
    filter_params = butter(5,0.03,btype='highpass')
    smoothed_data = filtfilt(filter_params[0],filter_params[1],raw_data)
    smoothed_data = gaussian_filter(smoothed_data,2.0)
    thresholded_data = smoothed_data > sigmas_smoothed[chname]*threshold_sig
    threshold_edges = np.convolve([1, -1], thresholded_data, mode='same')
    thresholded_edge_indices = np.where(threshold_edges==1)[0]
    
    return thresholded_edge_indices

#find peak in one event and display
#evt_num = 3
evt = df05.iloc[evt_num]

plt.rcParams['figure.figsize'] = [12,12]
plt.rcParams.update({'font.size': 18})

chlabels = []
yvalues = []
num_channels = 0

for i in range(len(evt['Channels'])):
    if 'SiPM' not in evt['ChannelTypes'][i]:
        continue
    else:
        num_channels += 1
    
    
    ch_name = analysis_config.GetChannelNameForSoftwareChannel(evt['Channels'][i])
    
    chlabels.append(ch_name)
    yvalues.append((num_channels-1)*500.)
    
    baseline = np.mean(evt['Data'][i][0:200]) 
    
    raw_data = evt['Data'][i]-baseline
    xvals = np.arange(len(raw_data))    
    plt.plot(xvals, raw_data + yvalues[num_channels-1],\
             label='Ch {}'.format(ch_name),color=color_cycle[num_channels-1])
    
    if '2-1'  not in ch_name and '2-3' not in ch_name and '2-4' not in ch_name:    
    #if '1-1' not in ch_name and '1-3' not in ch_name:
        pulse_idxs = PulseFinder(raw_data,ch_name,threshold_sig=4)
        for idx in pulse_idxs:
            plt.plot(xvals[idx], np.max(raw_data[idx:idx+20]) + yvalues[num_channels-1],'ok')
    
    
plt.xlim(-100.,12000.)    
plt.legend(fontsize=18,loc='upper right',handlelength=0.75)
plt.xlabel('Time (samples)')
plt.ylabel('Voltage ADC counts (500ADC offset)')
plt.title('Raw waveforms')
plt.savefig('./peak_display_one_event.png')
plt.close()

#SiPM gain calibration process 
pulse_heights = dict()

for label in chlabels:
    if '2-1'  not in label and '2-3' not in label and '2-4' not in label:    
    #if '1-1' not in label and '1-3' not in label:
        pulse_heights[label] = []
    


for index, evt in df05.iterrows():
    
    if index % 100 == 0: 
        print('Running event {}'.format(index))
    
    for i in range(len(evt['Channels'])):
        if 'SiPM' not in evt['ChannelTypes'][i]:
            continue

        ch_name = analysis_config.GetChannelNameForSoftwareChannel(evt['Channels'][i])

        if '2-1'  not in ch_name and '2-3' not in ch_name and '2-4' not in ch_name:    
        #if '1-1' in ch_name or '1-3' in ch_name:

            baseline = np.mean(evt['Data'][i][0:200]) 
            raw_data = evt['Data'][i]-baseline
        
            evtpulseidxs = PulseFinder(raw_data,ch_name,threshold_sig=3)#by yasheng
        
            evtpulseheights = []
            for idx in evtpulseidxs:
                evtpulseheights.append( np.max(raw_data[idx:idx+20]) )
        
            pulse_heights[ch_name].extend( evtpulseheights )
 
            #print(evtpulseheights)

for channel, data in pulse_heights.items():
    pulse_heights[channel] = np.array(data)


#Draw Calibration results and fix it with multigauss function
def SPEModel(x, An, mun, sign, Aspe, muspe, sigspe, Adpe, Atpe):
    noise = An * np.exp( -(x-mun)**2/(2*sign**2))
    signal = Aspe * np.exp( -(x-muspe)**2/(2*sigspe**2)) + \
        Adpe * np.exp( -(x-2*muspe)**2/(2*sigspe**2)) + \
        Atpe * np.exp( -(x-3*muspe)**2/(2*sigspe**2))
    return noise + signal

plt.rcParams['figure.figsize'] = [12,12]
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(ncols=2,nrows=3,sharex=True,sharey=True,gridspec_kw={'hspace':0.,'wspace':0.})



pulseheight_hists = dict()
for channel, data in pulse_heights.items():
    pulseheight_hists[channel] = hl.hist(pulse_heights[channel],bins=np.linspace(0.,1500.,150))

i = 0
for channel, hist in pulseheight_hists.items():
    
    row = (i) % 3
    col = int(np.floor( i/3. ))

    bin_centers = (hist.bins[0][1:] + hist.bins[0][0:-1])/2.
    bin_vals = hist.values
    bin_err = np.sqrt(hist.values)
    
    fitmin = 70.
    fitmax = 500.
    speguess = 150.
    sigguess = 30.
    aguess = 600.
    if '2-2' in channel:
        continue
    if '1-4' in channel:
        fitmax=600.
        fitmin=70.
        speguess=210.
        sigguess=60.
        aguess =1000.
    
    fitmask = (bin_centers<fitmax) & (bin_centers > fitmin)
    
    #p,pcov = opt.curve_fit(SPEModel,bin_centers[fitmask],bin_vals[fitmask],sigma=bin_err[fitmask],\
    #                      p0=(20000.,70.,10.,aguess,speguess,sigguess,200.,50.))
    
    #print('SPE mean in channel {}: {:4.4}'.format(channel,p[4]))
    
    hl.plot1d( ax[row,col], hist, color=color_cycle[i], label='data'.format(channel) )
    
    #xfit = np.linspace(0.,600.,300)
    #yfit = SPEModel(xfit,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7])
    #ax[row,col].plot(xfit,yfit,'--',color=color_cycle[i],\
                    #linewidth=1,label='Model fit\nSPE={:4.4} ADC'.format(p[4]))
    #yspe = Gaussian(xfit,p[3],p[4],p[5])
    #ydpe = Gaussian(xfit,p[6],2*p[4],p[5])
    #ytpe = Gaussian(xfit,p[7],3*p[4],p[5])

    #ax[row,col].plot(xfit,yspe,'--',color=color_cycle[i],linewidth=1)
    #ax[row,col].plot(xfit,ydpe,'--',color=color_cycle[i],linewidth=1)
    #ax[row,col].plot(xfit,ytpe,'--',color=color_cycle[i],linewidth=1)
    
    ax[row,col].text(1230.,800.,'{}'.format(channel),fontsize=20,\
tter                     bbox={'facecolor':(1.,1.,1.), 'edgecolor':(0.,0.,0.), 'alpha':1.})
    
    ax[row,col].legend(fontsize=15,edgecolor=(1.,1.,1.),framealpha=0.)
    ax[row,col].set_xlabel('Pulse height (ADC)')
    ax[row,col].set_yscale('log')
    #ax[row,col].set_ylim(5.,np.max(bin_vals)*1.1)
    ax[row,col].set_xlim(0.,1500.)

    i+=1
ax[1,0].set_ylabel('Counts')

plt.savefig('./spe_calibration_triggering_channels.png')#,dpi=200)#,bbox_inches='tight')












print('everything is ok!')
