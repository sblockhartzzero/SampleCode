function [y, waveform] = gen_Schroeder_Phase_Waveform( low_Freq, freq_Increment, num_Harmonics, schroeder_Sense, one_Window, T )
% function y = gen_Schroeder_Phase_Waveform( low_Freq, freq_Increment, num_Harmonics, glide_Direction, one_Window, T )
% 
% This function generates one pulse (aka burst) of a Schroder-phase
% waveform, given the following input parameters:
%       -low_Freq: the lowest frequency in Hz
%       -freq_Increment: the increment from one frequency component to
%       another aka the fundamental frequency
%       -num_Harmonics: the number of harmonics (not including the
%       component at low_Freq)
%       -schroeder_Sense is either "+" or "-", indicating +Schroeder or -Schroeder 
%       -one_Window: one window (column vector) which represents the
%       envelope.  It's max amplitude should be 1.
%       -T = 1/Fs i.e. the sample interval

% Determine the length of the window.
num_Points = length(one_Window);

% Generate the waveforms
% i.e. calculate the low_Freq component and all the harmonics in a waveform matrix.  
% (Later, we'll add them all up.)
waveform = zeros(num_Harmonics+1,num_Points);
for m = 1:num_Harmonics+1
    % Calculate the frequency of this component (harmonic)
    this_Freq = low_Freq + (m-1)*freq_Increment;
    % Calculate the phase offset of this component (according to the
    % Schroeder phase rule
    if strcmp(schroeder_Sense,'+')
        % This is a +Schroeder waveform, where m=2 is the 1st harmonic
        phase_Offset = pi*(m-1)*(m-2)/num_Harmonics;
    else
        phase_Offset = -pi*(m-1)*(m-2)/num_Harmonics;
    end    
    %phase_Offset = 0;
    for n = 1:num_Points
        waveform(m,n) = sin((2*pi*this_Freq*(n-1)*T) + phase_Offset );
    end
end

schroeder_Phase_Waveform = zeros(1,num_Points);
% Add up the components
for m = 1:num_Harmonics+1
    schroeder_Phase_Waveform = schroeder_Phase_Waveform + waveform(m,:);
end    

% Apply window by forming the dot product of the window and the carrier.
schroeder_Phase_Waveform_Windowed = one_Window' .* schroeder_Phase_Waveform;

% Return
y = schroeder_Phase_Waveform_Windowed;

