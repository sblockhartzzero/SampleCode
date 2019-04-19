function [mixture_Reconstructed, binary_Mask_Pattern, gamma_Signal, gamma_Masker, cf] = ideal_Binary_Mask(signal, masker, Fs, SNR_Threshold, CFG)

% Input Args
%   signal:                 Row-vector of signal
%   masker:                 Row-vector of masker segment (of same length as signal)
%   Fs:                     Sample rate (in Hz) of signal and masker
%   SNR_Threshold:          Threshold in dB e.g. -6 or -Inf
%   CFG:                    A structure containing cfg parameters
%                               -num_Filters (# of filters in the filter bank)
%                               -window_Duration_ms (duration of 50% overlapping windows in t-f plane, in ms)
%                               -min_Center_Frequency (lowest center frequency, in Hz)
%                               -max_Center_Frequency (highest center frequency, in Hz)
%
% Outputs
%   mixture_Reconstructed:  Time-domain representation after ideal binary mask is applied to the signal and masker
%   binary_Mask_Pattern:    Array of dim #overlapping windows x #filters, showing pattern of ideal binary mask
%   gamma_Signal:           Array of dim #filters x #samples, showing signal after gammatone filter bank is applied
%   gamma_Masker:           Array of dim #filters x #samples, showing masker after gammatone filter bank is applied
%   cf:                     Array of center frequencies used for gammatone filters
%
% This program generates the time domain representation of the ideal binary mask.  
% The caller must first make sure that the signal and masker are of the same length and same sample rate, Fs.
%
% This function processes both signal and masker through a gammatone filter bank. 
% Within each band, the filtered signal is compared to the filtered masker
% in 50% overlapping segments.  (These segments are shaped by a
% COLA-compliant window.)  If the SNR in the window is greater than the
% threshold, that window is turned on.  
%
% To get back to the time-domain, we just add up the filtered signal + noise across the N filters to get the
% "binary masker".  (Reverse filtering is not needed, as the gammatone_c filter does not introduce phase distortion.  
% See the utility gammatone_Driver.m to verify this.)
%
% Dependencies:
%   -This program calls gammatone_c, an implementation of a 4th-order
%   gammatone filter by Ning Ma of the University of Sheffield.


%%
%
% PREP SIGNAL and MASKER
%

% Inspect signal
signal_max_ampl = max(abs(signal));
signal_Length = length(signal);

% Parms
num_Filters = CFG.num_Filters;                              % Number of gammatone filters in the filter bank
window_Duration_ms = CFG.window_Duration_ms;                % Duration in ms of the sliding window
overlap = 0.5;                                              % Fraction of overlap between consecutive sliding windows
std_Threshold = 10^(SNR_Threshold/20);                      % Threshold for binary mask, compared to std

%
% Create overlapping windows
%
% If window length is odd, make it even.  This will help to make sure COLA
% (Constant Overlap Add) reconstruction doesn't have any problems.
window_Length = floor(Fs*window_Duration_ms/1000);          % Number of samples in the sliding window
if mod(window_Length,2) > 0 
    window_Length = window_Length + 1;
end
% Determine number of overlapping windows
% window_Start and window_End are used to figure out the start and end of
% the num_Windows overlapping windows.  If necessary, the time series is
% padded with zeros so that we have an integer number of overlapping
% windows.  Although the code below strictly supports any amount of
% overlap, to maintain COLA rules, we will only support 50% overlap in
% practice.
window_Start(1) = 1;
window_End(1) = window_Start(1) + window_Length - 1;
max_Val = window_End(1);
num_Windows = 1;
while max_Val < signal_Length
    % Next window
    num_Windows = num_Windows + 1;
    % Next window starting position
    window_Start(num_Windows) = window_End(num_Windows-1) - floor(overlap*window_Length) + 1;
    % Next window ending position
    window_End(num_Windows) = window_Start(num_Windows) + window_Length - 1;
    % Reset max_Val
    max_Val = window_End(num_Windows);
end
% Pad the signal for the last overlapping window, so we have an
% integral number of overlapping windows
length_Of_Padding = window_End(num_Windows) - signal_Length;
padding = zeros(1,length_Of_Padding);
signal_Padded = [signal,padding];
signal_Padded_Length = length(signal_Padded);

% Pad the masker as well
masker_Padded = [masker,padding];



%%
%
% FILTER
%
% Center freqs of filter bank
%cf = logspace(log10(CFG.min_Center_Frequency),log10(CFG.max_Center_Frequency),num_Filters);

cf = zeros(1,num_Filters);
cf(1) = CFG.min_Center_Frequency;
% The following formula comes from  a) cf1 + 0.5*BW1 = cf2 - 0.5*BW2 ...and
%                                   b) BW1 = 1.019*ERB1 (for gammatone filter) ...and
%                                   c) ERB1 = 24.7 + 0.108*cf1
for k = 2:num_Filters 
    cf(k) = ( (1 + 0.5*1.019*0.108)*cf(k-1) + (1.019*24.7) )/(1 - 0.5*1.019*0.108);
end


% Initialize filter outputs
gamma_Masker = zeros(num_Filters,signal_Padded_Length);
gamma_Signal = zeros(num_Filters,signal_Padded_Length);

% Call gammatone_c.
for k = 1:num_Filters
    gamma_Signal(k,:) = gammatone_c(signal_Padded, Fs, cf(k));
    gamma_Masker(k,:) = gammatone_c(masker_Padded, Fs, cf(k));
end

%%
%
% APPLY BINARY MASK
%
% Initialize binary_Mask_TF (filtered signal that is over threshold--still in t-f domain)
binary_Mask_TF = zeros(num_Filters,signal_Padded_Length);
binary_Mask_TF_flipped = zeros(num_Filters,signal_Padded_Length);
binary_Mask_TF_flipped_refiltered = zeros(num_Filters,signal_Padded_Length);

% Create window that is COLA-compliant
w = hann(window_Length, 'periodic');

% For each filtered time series, identify which time slices to keep
binary_Mask_Pattern = zeros(num_Windows,num_Filters);   
for k = 1:num_Filters
    for j = 1:num_Windows
        std_Sig = std( w'.*gamma_Signal(k,window_Start(j):window_End(j)) );
        std_Msk = std( w'.*gamma_Masker(k,window_Start(j):window_End(j)) );
        if ( std_Sig/std_Msk) >= std_Threshold
            % The SNR threshold has been statisfied, so retain this t-f window
            binary_Mask_Pattern(j,k) = 1.0;
        end
    end
end

% Build mask in t-f (binary_Mask_TF)
for k = 1:num_Filters
    for j = 1:num_Windows
        if binary_Mask_Pattern(j,k) == 1.0
            % Add the overlapping window to binary_Mask_TF
            m = 1;
            for n = window_Start(j):window_End(j)
                binary_Mask_TF(k,n) = binary_Mask_TF(k,n) + w(m)*(gamma_Signal(k,n)+gamma_Masker(k,n));
                m = m + 1;
            end
        end
    end
end

%%
% 
% RECONSTRUCT TO TIME DOMAIN and SCALE
%
% Reconstruct in time domain
%{
for k = 1:num_Filters
    binary_Mask_TF_flipped(k,:) = fliplr(binary_Mask_TF(k,:));
    binary_Mask_TF_flipped_refiltered(k,:) = gammatone_c(binary_Mask_TF_flipped(k,:), Fs, cf(k)); 
end
mixture_Reconstructed_array = fliplr(binary_Mask_TF_flipped_refiltered);
mixture_Reconstructed = sum(mixture_Reconstructed_array);
%}
% Instead of reverse filtering, try just adding up filtered time series 
mixture_Reconstructed_array = binary_Mask_TF;
mixture_Reconstructed = sum(mixture_Reconstructed_array);

% Scale recon--THIS IS A BUG!!
%mixture_Reconstructed_max_ampl = max(abs(mixture_Reconstructed));
%signal_recon_scaled = (signal_max_ampl/mixture_Reconstructed_max_ampl)*mixture_Reconstructed;
% INSTEAD, DO NOT SCALE THE OUTPUT!!!




