function [Fs, mean_Array, max_Array, var_Array, datenum_Array, event_Whistle] = whistle_Driver(input_File_Fullpath, experiment, signal_Type,start_Time_Datenum, WHISTLE_CFG)

% This program reads the individual wav file specified by input_File_Fullpath and calls Silbido to get whistle events for this file. Like the qc_Driver program, 
% it also gathers statistics (mean, variance) per 50%-overlapping window, where the window size is specified in the WHISTLE_CFG structure. It does
% this in order to estimate the IN-BAND noise. So, before measuring the variance per window, the entire time series is high-pass filtered, given
% the cutoff frequency in WHISTLE_CFG.fco_HP

%{
INPUTS:
    -input_File_Fullpath:                           Full path of the sound file, ending with a filesep
    -experiment:                                    e.g. CSI1 to indicate first deployment of CSI dataset
    -signal_Type:                                   Currently either wav or aif
    -start_Time_Datenum:                            Datenum for the start of the sound file
    -WHISTLE_CFG
        -WHISTLE_CFG.channel:                       If the sound file (e.g. wav) stores multiple channels, you can specify which channel to analyze
        -WHISTLE_CFG.analysis_Window_Duration_Secs: Duration (in secs) of the window for QC analysis
        -WHISTLE_CFG.plots:                         Boolean. If true, generate plot(s) per set
        -WHISTLE_CFG.fco_Hz_HP                      Cut-off frequency (in Hz) for high-pass filter (used only to estimate the IN-BAND noise)


OUTPUTS
    -Fs:                                            The sample rate of this sound file
    -mean_Array:                                    Row vector of the mean values (per window) of the signal, of dimension 1 x #windows
    -max_Array                                      Row vector of the max (abs) values (per window) of the signal, of dimension 1 x #windows
    -var_Array:                                     Row vector of the variances (per window) of the signal, of dimension 1 x #windows
    -datenum_Array:                                 Row vector of the midpoint of the window (in datenum format), of dimension 1 x #windows 
    -event_Whistle:                                 Struct with attributes start_Datenum, stop_Datenum, as well as statistics on the frequencies
                                                        and snr
%}

% 20170130:     Changed SNR thrshold to 10 (from 7), both here in call to
% dtTonalsTracking and in C:\Users\NewFolderSamsung\Share\Silbido\matlab\spectral_analysis\dtThresh.m

input_File_Fullpath


%% INIT RETURN VARS
% event_Whistle
event_Whistle = struct;


%% PARMS
secs_Per_Day = 3600*24;


%% GET SIGNAL, SAMPLE RATE
switch signal_Type
    case 'wav'
        [input_Data,Fs] = wavread(input_File_Fullpath);                     % signal is a column vector
        [num_Samples,num_Channels] = size(input_Data);
        % Set the signal to the specified input channel
        if WHISTLE_CFG.channel > num_Channels
            error_String = sprintf('%s %s %s %s','The signal has only',int2str(num_Channels),'channels, file = ',input_File_Fullpath);
            error(error_String);
        else
            signal = input_Data(:,WHISTLE_CFG.channel);
        end
    case 'aif' 
        [input_Data,Fs,wmode,fidx]=readaif(input_File_Fullpath);            % signal is a column vector
        [num_Samples,num_Channels] = size(input_Data);
        % Set the signal to the specified input channel
        if WHISTLE_CFG.channel > num_Channels
            error_String = sprintf('%s %s %s %s','The signal has only',int2str(num_Channels),'channels, file = ',input_File_Fullpath);
            error(error_String);
        else
            signal = input_Data(:,WHISTLE_CFG.channel);
        end
    otherwise
        error('Unknown signal type');
end

%% HP FILTER
% Filter first
Fc = WHISTLE_CFG.fco_Hz_HP;
Fc_Normalized = Fc/(Fs/2);
[b,a]=butter(4,Fc_Normalized,'high');
signal_Filtered = filtfilt(b, a, signal);


%% GET STATS PER WINDOW
% Get dimensions
window_Length_Samples = floor(WHISTLE_CFG.analysis_Window_Duration_Secs*Fs) - 1;
if num_Samples < window_Length_Samples
    % Handle special case
    num_Windows = 1;
else
    % Calc number of 50% overlapping windows
    num_Windows = 2*(floor(num_Samples/window_Length_Samples)) - 1;
end
% Init arrays
mean_Array = zeros(1,num_Windows);
max_Array = zeros(1,num_Windows);
datenum_Array = zeros(1,num_Windows);
var_Array = zeros(1,num_Windows);
% Stats per window
for window_Num = 1:num_Windows
    % Get start, stop, and mid of this window
    if num_Windows == 1
        start_Sample = 1;
        stop_Sample = num_Samples;
    else
        start_Sample = 1 + floor((window_Num - 1)*window_Length_Samples/2);
        stop_Sample = start_Sample + window_Length_Samples;
    end
    mid_Sample = floor((stop_Sample + start_Sample)/2);
    % Time for mid-point of the window
    time_Window_Secs = mid_Sample/Fs;
    datenum_Array(window_Num) = start_Time_Datenum + (time_Window_Secs/secs_Per_Day);
    % Stats per window
    mean_Array(window_Num) = mean(signal_Filtered(start_Sample:stop_Sample));
    max_Array(window_Num) = max(abs(signal_Filtered(start_Sample:stop_Sample)));
    if max_Array(window_Num) == 0
        display_String = sprintf('%s %s %s %s','Dropout in file',input_File_Fullpath, 'for window',int2str(window_Num));
        disp(display_String);
    end
    var_Array(window_Num) = var(signal_Filtered(start_Sample:stop_Sample));
end



%% EVENTS
% For Whistles, we need to call silbido first. Caller must addpath to
% C:\Users\NewFolderSamsung\Share\Silbido and call silbido_init
% Try saving filtered signal to wav file before calling Silbido (as I'm not
% sure Silbido is filtering first).
wavwrite(signal_Filtered,Fs,'signal_Filtered.wav');
% Call silbido
%detections = dtTonalsTracking(input_File_Fullpath,0,Inf,'Threshold',5);
try
    detections = dtTonalsTracking('signal_Filtered.wav',0,Inf,  'Threshold',10);
%                                                            'Noise','none',...
%                                                            'RemoveTransients', true);
catch err
    warn_String = strcat('WARN: Call to Silbido failed, so skipping',input_File_Fullpath);
    disp(warn_String);
    event_Whistle.start_Datenum = [];
    event_Whistle.stop_Datenum = [];
    event_Whistle.f_Hz_Min = [];
    event_Whistle.f_Hz_Max = [];
    event_Whistle.f_Hz_Median = [];
    event_Whistle.snr_Power_dB_Median = [];
    return
end
% Convert java.util.LinkedList to array of java.lang.Object
detections_Array = detections.toArray();
num_Events = length(detections_Array);
% For each event, get the t,f associated with the min,max frequency. Also,
% determine the median frequency for the event
% Init
event_Whistle.start_Datenum = nan(1,num_Events);
event_Whistle.stop_Datenum = nan(1,num_Events);
event_Whistle.f_Hz_Min = nan(1,num_Events);
event_Whistle.f_Hz_Max = nan(1,num_Events);
event_Whistle.f_Hz_Median = nan(1,num_Events);
event_Whistle.snr_Power_dB_Median = nan(1,num_Events);
for event_Num = 1:num_Events
    t_Secs = detections_Array(event_Num).get_time();                % Method of tonals
    f_Hz = detections_Array(event_Num).get_freq();                  % Method of tonals
    snr_Power_dB = detections_Array(event_Num).get_snr();
    % Get median freq (in Hz)
    f_Hz_Median = median(f_Hz);
    snr_Power_dB_Median = median(snr_Power_dB);
    % Get max, min freqs (In Hz)
    [f_Hz_Min, ~] = min(f_Hz);
    [f_Hz_Max, ~] = max(f_Hz);   
    % Stuff into event struct
    event_Whistle.start_Datenum(event_Num) = start_Time_Datenum + (t_Secs(1)/secs_Per_Day);
    event_Whistle.stop_Datenum(event_Num) =  start_Time_Datenum + (t_Secs(end)/secs_Per_Day);
    event_Whistle.f_Hz_Min(event_Num) = f_Hz_Min;
    event_Whistle.f_Hz_Max(event_Num) = f_Hz_Max;
    event_Whistle.f_Hz_Median(event_Num) = f_Hz_Median;
    event_Whistle.snr_Power_dB_Median(event_Num) = snr_Power_dB_Median;
end



%% Write events to a file (as Raven Selections)
% Clipping
if num_Events > 0
    event_Type = 'WHISTLE';
    % Specify dir where selections file will be stored
    switch experiment
        case 'CSI1'
            selections_Dir = 'C:\Users\slockhar\Projects\CSI\Output\CSI1\Selections\Whistles\';
        case 'CSI2'
            selections_Dir = 'C:\Users\slockhar\Projects\CSI\Output\CSI2\Selections\Whistles\';
        otherwise
            error('Unknown experiment');
    end  
    % Need wav filename sans extension
    filesep_Loc = strfind(input_File_Fullpath,filesep);
    extension = strcat('.',signal_Type);
    ext_Loc = strfind(input_File_Fullpath,extension);
    start_Loc = filesep_Loc(end) + 1;
    end_Loc = ext_Loc - 1;
    wavFilename_Sans_Extension = input_File_Fullpath(start_Loc:end_Loc);
    save_Selection_File( event_Type, selections_Dir, wavFilename_Sans_Extension, start_Time_Datenum, Fs, event_Whistle );
end





