"""
convert all CSV files (e.g. ID01_epoch1.csv.xz) to the same sampling frequency

John M. O'Toole, CergenX
Started: 30-10-2022
last update: Time-stamp: <2022-11-09 18:29:56 (otoolej)>
"""
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import signal as sp
import matplotlib.pyplot as plt
import pathlib
from multiprocessing import Pool


def resample_eeg_files(eeg_dir, out_dir, fs_new=64, lpf_cutoff=30, proc_parallel=False):
    """ Convert all EEG files in .csv.xz (or .csv) format

    Parameters
    ----------
    eeg_dir: string
        directory contain the CSV EEG files
    out_dir: string
        directory to write the new files
    fs_new: int
        new sampling frequency (default 64) 
    lpf_cutoff: float
        low-pass filter cut-off (default 30)
    proc_parallel: bool
        speed-up processing by doing in parallel (if memory allows)
    """
    if proc_parallel:
        all_files = pathlib.Path(eeg_dir).glob('*.[csv.xz csv]*')
        all_args = [(fname, lpf_cutoff, fs_new, out_dir) for fname in all_files]

        with Pool() as pool:
            pool.starmap(convert_file, all_args)

    else:
        for fname in pathlib.Path(eeg_dir).glob('*.[csv.xz csv]*'):
            convert_file(fname, lpf_cutoff, fs_new, out_dir)

    

def convert_file(fname, fc, fs_new, out_dir):
    """read .csv file, re-sample, and save to new .csv file

    Parameters
    ----------
    fname: string
        file name of .csv file or .csv.xz (full path to file)
    fc: float
        cut-off frequency for low-pass anti-aliasing filter
    fs_new: int
        sampling frequency
    out_dir: string
        name of directory to save new file to
    """
    file_name = pathlib.Path(fname).name

    eeg_df = pd.read_csv(fname)
    
    ttime = eeg_df['time'].values
    fs = round(1 / (ttime[1] - ttime[0]))

    ch_names = eeg_df.columns[1:]

    print(f'\t|{file_name}\t| {fs}\t|')

    if fs != fs_new:
        x_resample_dt = dict.fromkeys(['time'] + list(ch_names))
       
        for ch in ch_names:
            x_filt = filt_lowpass(eeg_df[ch].values, fs, fc, 'FIR')

            x_resample_dt[ch] = sp.resample_poly(x_filt, fs_new, fs)

        ttime_new = np.arange(0, int(len(x_filt) / fs), (1 / fs_new))
        x_resample_dt['time'] = ttime_new
    
        eeg_new_df = pd.DataFrame(x_resample_dt)
            
    else:
        eeg_new_df = eeg_df

    eeg_new_df.to_csv(pathlib.Path(out_dir) / file_name, index=False, compression='infer')


    
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
        l_filt = 21        
        sos = sp.cheby2(l_filt, Wn=fc, rs=100, btype='lowpass', fs=fs, output='sos')    

        # zero-phase filter:
        xmean = np.nanmean(x)        
        y = sp.sosfiltfilt(sos, x - xmean)
        y += xmean

        if db_plot:
            plot_freqz(sos, Fs=fs, Nfreq=2048, fig_num=1, ba_or_sos='sos')
        
    else:
        raise Exception("Unknown filter type; should be either FIR or IIR")
        
    return y
    



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



def _add_trailing_slash(ddir):
    if not ddir.endswith('/'):
        return ddir + '/'
    return ddir
    
    
def parsing_cli_args(args=[]):
    """ parse command line arguments """
    parser = argparse.ArgumentParser(description="Downsample EEG files")

    parser.add_argument('--eeg_dir', type=str, required=True, 
                        help="directory holding the compressed CSV files")
    parser.add_argument('--out_dir', type=str, required=True, 
                        help="directory to write downsampled files to")
    parser.add_argument('--fs_new', type=float, default=64.0,
                        help="new sampling frequency (in Hz)")
    parser.add_argument('--lpf_cutoff', type=float, default=30.0,
                        help="cut-off frequency for anti-aliasing filter (in Hz)")
    parser.add_argument('--proc_parallel', type=bool, default=False,
                        help="process in parallel?")
    args = parser.parse_args()

    return args



if __name__ == "__main__":
    # 1. first, read command line arguments (if any)
    args = parsing_cli_args(sys.argv[1:])

    from pprint import pprint
    pprint(vars(args))

    resample_eeg_files(**vars(args))
    
    
    


    


    
