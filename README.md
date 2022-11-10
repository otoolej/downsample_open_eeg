# Downsampling routines for an open-access neonatal EEG dataset

Python and Matlab code to downsample neonatal EEG files stored at the Zenodo repository [DOI:
10.5281/zenodo.6587973](https://doi.org/10.5281/zenodo.6587973).

The EEG is recorded from newborn infants in the neonatal intensive care unit using 2
different EEG machines, the NicoletOne ICU Monitor (Natus, Middletion, WI, USA) and the
Neurofax EEG-1200 (Nihon Khoden, Tokyo, Japan). EEG was sampled at a rate of 256 Hz for
the NicoletOne and 200 Hz for the Neurofax. The open dataset is described in reference
[1](#references).

The EEG data in [DOI:10.5281/zenodo.6587973](https://doi.org/10.5281/zenodo.6587973) is
stored in 2 formats: European Data Format (EDF) and a plain-text comma-separated values
(CSV) format. The EDF version can be used in most EEG reviewing software, for example
using [EDFbrowser](https://www.teuniz.net/edfbrowser/). The CSV version is a convient
format for computer processing and the code presented here reads this format. All CSV
files are in a compressed format using the [XZ](https://tukaani.org/xz/) compression format.



Most quantitative analysis and machine learning analysis will downsample the data to a
lower frequency, such as 64 Hz. This is because very little EEG signal power will be
present at high (>30 Hz) frequencies in neonatal EEG and therefore downsampling can
greatly reduce the required computational resources without loosing significant signal
information.

We present downsampling routines to process the EEG files using either Matlab or Python

- [Python code](#python-code)
- [Matlab code](#matlab-code)


---

## Python code
Tested with Python 3.9 and 3.10. See [requirements.txt](python/requirements.txt) file for
a list of modules and the versions used for testing. 

In a shell, run:
```sh
$ cd python_code/
$ python3 -m resample_xz_files --eeg_dir=<LOCATION_OF_EEG_FILES> --out_dir=<EMPTY_DIRECTORY_FOR_NEW_FILES>
```

All available command-line parameters:

| parameter       | description                                            | default |
|-----------------|--------------------------------------------------------|---------|
| `eeg_dir`       | directory of the compressed CSV files; REQUIRED        |         |
| `out_dir`       | empty directory for new compressed CSV files; REQUIRED |         |
| `fs_new`        | new sampling frequency in Hz                           | 64      |
| `lpf_cutoff`    | cut-off frequency for anti-aliasing filter             | 30      |
| `proc_parallel` | do the processing in parallel?                         | False   |

Can also run directly within python interpreter/notebook by calling the
`resample_eeg_files` function within the file.

---

## Matlab code

Tested with Matlab R2020a. Requires the Signal Processing Toolbox.

All EEG files musts be decompressed first before runing these Matlab routines, as at time
of writing, Matlab lacks support to read files compressed with
[XZ](https://tukaani.org/xz/).

Change directory to `matlab_code` and within Matlab run as:

```matlab
>> resample_eeg_files(<LOCATION_OF_EEG_FILES>, <EMPTY_DIRECTORY_FOR_NEW_FILES>, 64, 30);
```

where 64 is the new sampling frequency and 30 is the low-pass filter cutoff frequency.

---


# References

1. JM O'Toole, SR Mathieson, SA Raurale, F Magarelli, WP Marnane G Lightbody & GB
   Boylan. Neonatal EEG graded for severity of background abnormalities in
   hypoxic-ischaemic encephalopathy. [arXiv preprint
   arXiv:2206.04420](https://arxiv.org/abs/2206.04420). 2022 Jun 9.
2. JM O'Toole, SR Mathieson, F Magarelli, WP Marnane, G Lightbody & GB
   Boylan. (2022). Neonatal EEG Graded for Severity of Background Abnormalities (1.0)
   [Data set]. Zenodo. [10.5281/zenodo.6587973](https://doi.org/10.5281/zenodo.6587973)
