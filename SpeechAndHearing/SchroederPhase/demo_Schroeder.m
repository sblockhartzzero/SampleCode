clear

% This script tests the generation of Schroeder-phase waveforms

% Prepare data for this run
Fs=5000;
signal_Duration_Secs = 0.70;


% Prep call to gen_Schroeder_Phase_Waveform
% Use values like in Wojtczak and Oxenham (2009)
 low_Freq = 400;
 freq_Increment = 1;
 num_Harmonics = 600;
 schroeder_Sense = '-';
 one_Window = ones(1,signal_Duration_Secs*Fs);
 T=1/Fs;
 
 % Call gen_Schroeder_Phase_Waveform
 [y, waveform] = gen_Schroeder_Phase_Waveform( low_Freq, freq_Increment, num_Harmonics, schroeder_Sense, one_Window', T );
 
 % Play
 %y_sch = [zeros(1,1000),y,zeros(1,1000)];
 y_sch = y/(1.1*max(abs(y)));
 sound(y_sch,Fs);
 
 % Plot
 figure; plot(y_sch);
 figure; specgram(y_sch,256,Fs);
 
 % Check
 %x = sum(waveform);
 %figure; specgram(x,256,Fs);
 
 % Save
 audiowrite('sampleSchroederPhase.wav',y_sch,Fs);
 
 
