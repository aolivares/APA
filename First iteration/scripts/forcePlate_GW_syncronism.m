% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------- GW and forceplate syncronism --------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Bötzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Verónica  Torres (2): 
%                   |_ vts24@correo.ugr.es 
%                 - Dr. Eng. Alberto Olivares (3): 
%                   |_ aolivares@ugr.es
%                 - Robin Weiss (4): 
%                   |_ mail@robinweiss.de
%
% * Affiliation: (1) Laboratory of Motion Analysis and Deep Brain 
%                    Stimulation, Department of Neurology, Klinikum 
%                    Grosshadern, Munich, Germany.
%                (2) Master in Telecommunication Engineering, University of 
%                    Granada, Granada, Spain, (student).
%                (3) Signal Processing and Biomedical Applications group,
%                    Department of Signal Theory, Telematics and
%                    Communications, University of Granada, Granada, Spain.
%                (4) Bachelor in Electrical Engineering, University of 
%                    Applied Sciences of Munster, Munster, Germany, 
%                    (student).
%
% * Last modification: 08/01/2015
% -------------------------------------------------------------------------
% INFORMATION: This file contain the routine to to obtain the differents 
% signals (differents repeats) from GW and synchronize them. After, these 
% signals are sinchronized with AP excursions from forceplate.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clear workspace and close all figures.
% -------------------------------------------------------------------------
clear all; close all; clc;

% -------------------------------------------------------------------------
% Import GaitWatch functions library.
% -------------------------------------------------------------------------

gw = gwLibrary;

% -------------------------------------------------------------------------
% Load data from GaitWatch and forceplate
% -------------------------------------------------------------------------

% Load the data from ES39_05 (forceplate) and GW_1605_150814_1154 (GW)
% We use pitch of rigth shank signal from GW and both feet AP from 
% foceplate to synchronize.

load ('forcePlate_GW_data_ES39\both_AP');
AP_forceplate=both_AP;
load ('forcePlate_GW_data_ES39\timeFP');
time_forceplate=timeFP;
load ('forcePlate_GW_data_ES39\pitch_GKF_right_shank');
pitch_GW=pitch_GKF_right_shank;
load ('forcePlate_GW_data_ES39\time');
time_GW=time;

% -------------------------------------------------------------------------
% Detect the first peak in the forcepate signal
% -------------------------------------------------------------------------

% Find all peaks greater than 100.
[peak_value, peak_location] = findpeaks(AP_forceplate,'minpeakheight',100);

% We use the first peak to synchronize.
peak_loc_FP=peak_location(1);

% -------------------------------------------------------------------------
% Separate the differents episodes and deffine the peak off of GW signal
% -------------------------------------------------------------------------

% We obtain the differents signals from GW in only one matriz and the index
% used to synchronism (peak off location).
[signal_separated, peak_loc_GW]=Separate_episodes(pitch_GW);

% Define the new time array for these signals. It's a approximation because
% it's really the longest interval selected previously and used to synchronize.
time_signals=time_GW(1:length(signal_separated(1,:)));

% -------------------------------------------------------------------------
% Plot GW signals
% -------------------------------------------------------------------------

figure (1)
plot(signal_separated');
title('Pitch of rigth shank');
xlabel('Samples');
ylabel('Pitch (deg)');
text(peak_loc_GW, -35, '\uparrow', 'VerticalAlignment', 'Bottom', ...
    'HorizontalAlignment','Center','FontSize',16,'col','r');

figure (2)

plot(time_signals,signal_separated');
title('Pitch of rigth shank');
xlabel('Time (s)');
ylabel('Pitch (deg)');

% -------------------------------------------------------------------------
% Interpolation of forceplate signal.
% -------------------------------------------------------------------------

% In this case, we are going to interpolate with the most duration signal.
% Set the vector with the points that we need to interpolate.
div=peak_loc_FP/peak_loc_GW;
xi1=1:div:peak_loc_FP;

% Compute a linear interpolation between the begining and the peak
% location.
AP_FP1=interp1(1:peak_loc_FP,AP_forceplate(1:peak_loc_FP),...
    xi1,'spline');

% Campute a linear interoplation between the peak location and the end.
div=(length(AP_forceplate)-peak_loc_FP)/(length(signal_separated(1,:))-peak_loc_GW);
xi2=peak_loc_FP:div:length(AP_forceplate);

AP_FP2=interp1(peak_loc_FP:length(AP_forceplate),...
    AP_forceplate(peak_loc_FP:length(AP_forceplate)),xi2,'spline');

% Join both interpolations

AP_FP_interpolated=[AP_FP1 AP_FP2];

% -------------------------------------------------------------------------
% Plot the first GW signal and forceplate signal
% -------------------------------------------------------------------------

figure (3)
subplot(3,1,1)
    plot(AP_FP_interpolated);
    title('Ant-Post COP excursions');
    xlabel('Samples');
    ylabel('AP-COP (mm)');
    
subplot(3,1,2)
    plot(signal_separated(1,:),'r');
    title('Pitch of rigth shank');
    xlabel('Samples');
    ylabel('Pitch (deg)');

subplot(3,1,3)
    plot(AP_FP_interpolated)
    hold on
    plot(signal_separated(1,:),'r')
    legend('AP excursions (forceplatform)','pitch rigth shank (GaitWatch)');
    title('Synchronized Signals (GW and forceplatform)');
    hold off




% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END 
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\