%-------------------------------------------------------------------------------
% resample_eeg_files: 
%
% Syntax: [] = resample_eeg_files(eeg_dir, out_dir)
%
% Inputs: 
%     eeg_dir, out_dir - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 30-10-2022
%
% last update: Time-stamp: <2022-10-30 21:47:55 (otoolej)>
%-------------------------------------------------------------------------------
function resample_eeg_files(eeg_dir, out_dir)

%---------------------------------------------------------------------
% PARAMETERS: set here
%---------------------------------------------------------------------
% low-pass cutoff for antialiasing filter:
fc = 95;
% new sampling frequency:
fs_new = 200;


%---------------------------------------------------------------------
% locate files and process each individually
%---------------------------------------------------------------------
lfiles = dir([eeg_dir '*.csv']);

if(isempty(lfiles))
    error('No .csv files found. Make sure folder is correct and files are uncompressed.');
    
else
    % iterate over all files and save resampled file to out directory:
    for n = 1:length(lfiles)
        fname = [lfiles(n).folder filesep lfiles(n).name];
        
        convert_file(fname, fc, fs_new, out_dir);
    end
end



function convert_file(fname, fc, fs_new, out_dir)
%---------------------------------------------------------------------
% read .csv file, re-sample, and save to new .csv file
%---------------------------------------------------------------------

[~, file_name, fext] = fileparts(fname);
if(strcmp(fext, '.xz'))
    error('Need to uncompress .xz files first before opening in Matlab.');
else

    eeg_tb = readtable(fname);
    
    ttime = eeg_tb.time;
    fs = round(1 / (ttime(2) - ttime(1)));
    
    ch_names = eeg_tb.Properties.VariableNames(2:end);
    
    
    fprintf('\t|%s\t| %d\t| \n', file_name, fs);

    
    if(fs ~= fs_new)
        x_resample = [];
        for n = 1:length(ch_names)
            x = eeg_tb{:, ch_names{n}};
            
            % apply a low-pass anti-aliasing filter:
            x_filt = filt_lowpass(x, fs, fc, 'FIR');

            % up or downsample:
            x_resample{n + 1} = resample(x_filt, fs_new, fs);
        end
        ttime_new = (0:(1 / fs_new):(length(x) / fs));
        x_resample{1} = ttime_new([1:end - 1])';
        
        eeg_new_tb = table(x_resample{:}, 'VariableNames', [{'time'}, ch_names]);
    else
        eeg_new_tb = eeg_tb;
    end

    
    % write to new file
    fname_out = [out_dir, file_name, '.csv'];
    writetable(eeg_new_tb, fname_out);
end




function y = filt_lowpass(x, fs, fc, filt_type)
%---------------------------------------------------------------------
% low pass filter
%---------------------------------------------------------------------
if(nargin < 4 || isempty(filt_type)), filt_type = 'FIR'; end


db_plot = false;

switch lower(filt_type)
  case 'fir'

    l_filt = 2000;
    w = blackmanharris(l_filt + 1);
    b = fir1(l_filt, fc / (fs / 2), w);

    xmean = nanmean(x);
    y = filtfilt(b, 1, x - xmean);
    y = y + xmean;
    
    if(db_plot)
        freqz(b, 1, 2048, fs);
    end


  case 'iir'
    %  if IIR 
    l_filt = 11;
    [z, p, k] = butter(l_filt, fc / (fs / 2), 'low');    
    
    [sos, g] = zp2sos(z, p, k);
    y = filtfilt(sos, g, x);

    if(db_plot)
        fvtool(sos, 'Fs', fs); 
    end

    
  otherwise 
    error("filt_type should be either FIR or IIR");
    
end





