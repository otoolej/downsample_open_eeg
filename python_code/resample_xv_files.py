"""
convert all CSV files (e.g. ID01_epoch1.csv.xz) to the same sampling frequency

John M. O'Toole, CergenX
Started: 30-10-2022
last update: Time-stamp: <2022-10-30 22:04:32 (otoolej)>
"""
import numpy as np
import pandas as pd
from scipy import signal as sp
import matplotlib.pyplot as plt


# -------------------------------------------------------------------
#  PARAMETERS, set here
# -------------------------------------------------------------------
FS_NEW = 200      # new sampling frequnecy for all files
LPF_CUTOFF = 95   # low-pass filter cut off



def filt_lowpass(x, fs, fc, filt_type='FIR', db_plot=False):
    """Low pass filter

    Parameters
    ----------
    x: ndarray
        input signal
    fs: integer
        sampling frequency
    fc: float
        low-pas cutoff frequency
    filt_type: string
        filter type, either 'FIR' or 'IIR' (default IIR)
    """
    if filt_type.lower() == 'fir':
        # design low-pass filter:
        l_filt = 2001
        b = sp.firwin(l_filt, fc, window='blackmanharris', pass_zero='lowpass', fs=fs)

        # zero-phase filter:
        xmean = np.nanmean(x)
        y = sp.filtfilt(b, 1, x - xmean, padlen=3 * (len(b) - 1))
        y += xmean

        if db_plot:
            plot_freqz(b, Fs=fs, Nfreq=2048, fig_num=1, ba_or_sos='ba')
        

    elif filt_type.lower() == 'iir':
        #  design high-pass Butterworth filter:
        l_filt = 11        
        sos = sp.butter(l_filt, fc, btype='lowpass', fs=fs, output='sos')    

        # zero-phase filter:
        xmean = np.nanmean(x)        
        y = sp.sosfiltfilt(sos, x - xmean)
        y += xmean

        if db_plot:
            plot_freqz(sos, Fs=fs, Nfreq=2048, fig_num=1, ba_or_sos='sos')
        
    else:
        raise Exception("Unknown filter type; should be either FIR or IIR")
        
        

    

def read_downsample_edf(fname, params=None, db_plot=False):
    """ 
    Parameters
    ----------
    fname: string
        file name of EDF file
    params: dataclasses
        parameters 
    db_plot: bool
        plot the EEG

    Returns
    -------
    eeg_down : ndarray
        downsampled EEG data as bipolar montage
    ch_refs : list
        channel names
    fs_new : int
        frequency of downsampled EEG
    imp_mask : ndarray
        1D mask to indicate impedance check (sampled at fs_new)
    """
    if params is None:
        params = ibi_parameters.ibiParams()
        

    # -------------------------------------------------------------------
    #  read in the EEG data (using the channels needed)
    # -------------------------------------------------------------------
    (eeg_data, fs, ch_refs, times, edf_raw, idx, start_date) = load_edf(fname, params, db_plot)
    n_channels, n_x = eeg_data.shape


    # -------------------------------------------------------------------
    #  filter each channel separately
    # -------------------------------------------------------------------
    # design the FIR filter:
    b = firwin(
        params.l_filt, params.fc, window=params.win_type, pass_zero="lowpass", fs=fs
    )
    
    for n in range(n_channels):
        # do the filtering:
        x = _do_lpf(eeg_data[n, :], b)

        db_plot_tmp = False
        if db_plot_tmp:
            plt.figure(2, clear=True)
            plt.plot(eeg_data[n, :])
            plt.plot(x)
            plt.show()
            plt.pause
        
        eeg_data[n, :] = x

    # -------------------------------------------------------------------
    #  downsample
    # -------------------------------------------------------------------
    x_down = _do_downsample(eeg_data[0, :], fs, params.fs_new)
    eeg_down = np.zeros((n_channels, len(x_down)))
    eeg_down[0, :] = x_down
    for n in range(1, n_channels):
        eeg_down[n, :] = _do_downsample(eeg_data[n, :], fs, params.fs_new)




    
    
def _do_downsample(x, fs, fs_new):
    """ do the downsampling for 1D signal """
    f_ratio = fs / fs_new
    
    if f_ratio == int(f_ratio):
        x_down = x[::int(f_ratio)]

    else:
        # if f_ratio is not an integer, then up sample then downsample to the correct
        # sampling rate:
        x_down = resample_poly(x, fs_new, fs, padtype='line')

    return x_down


def plot_freqz(b, a=1, Fs=1, Nfreq=None, fig_num=1, ba_or_sos='ba'):
    """plot magnitude and phase of filter

    Parameters
    ----------
    b: ndarray
        filter coefficients
    a: ndarray
        filter coefficients
    Fs: scalar
        sampling frequency
    Nfreq: scalar
        frequency-domain sampling frequency
    fig_num: scalar
        figure number
    """
    if ba_or_sos == 'ba':
        w, h = sp.freqz(b, a, worN=Nfreq)
    else:
        w, h = sp.sosfreqz(b, worN=Nfreq)

    ffr = (w / np.pi) * (Fs / 2)
    fig, ax1 = plt.subplots(num=fig_num, clear=True)
    plt.title('filter response')
    plt.plot(ffr, 20 * np.log10(abs(h) + np.finfo(float).eps), 'b')
    plt.ylabel('Amplitude [dB]', color='b')
    plt.xlabel('Frequency [Hz]')

    plt.twinx(ax1)
    angles = np.unwrap(np.angle(h))
    plt.plot(ffr, angles, 'g')
    plt.ylabel('Angle (radians)', color='g')
    plt.grid()
    plt.axis('tight')

    plt.pause(0.0001)
    plt.show()
