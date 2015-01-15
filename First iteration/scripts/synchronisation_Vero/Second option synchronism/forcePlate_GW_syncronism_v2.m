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

load ('forcePlate_GW_data_ES39\both_force');
force_forceplate=both_force;
load ('forcePlate_GW_data_ES39\timeFP');
time_forceplate=timeFP;
load ('forcePlate_GW_data_ES39\a_Z_right_shank_1_C');
acc_GW=a_Z_right_shank_1_C;
load ('forcePlate_GW_data_ES39\time');
time_GW=time;

% -------------------------------------------------------------------------
% When pacient set on foot over forceplate,we have a value in the signal,
% so the first sample is the beginning, an we must match with the first
% peak of the GW signals.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Separate the differents episodes and deffine the peak off of GW signal
% -------------------------------------------------------------------------

% We obtain the differents signals from GW in only one matriz and the index
% used to synchronism (peak off location).
[signal_separated, peak_loc_GW]=Separate_episodes_v2(acc_GW);

% Define the new time array for these signals. It's a approximation because
% it's really the longest interval selected previously and used to synchronize.
time_signals=time_GW(1:length(signal_separated(1,:)));

% -------------------------------------------------------------------------
% Plot GW signals
% -------------------------------------------------------------------------

figure (1)
subplot(2, 1, 1)
plot(force_forceplate);
subplot(2, 1, 2)
plot(acc_GW)

figure (2)
plot(signal_separated');
title('Acceleration of rigth shank');
xlabel('Samples');
ylabel('Acc');
text(peak_loc_GW, 0.2, '\uparrow', 'VerticalAlignment', 'Bottom', ...
    'HorizontalAlignment','Center','FontSize',16,'col','r');

figure (3)
subplot(2,1,1)
    plot(force_forceplate);
    title('Force of both feet');
    xlabel('Samples');
    ylabel('Force (N)');
    
subplot(2,1,2)
    plot(signal_separated(1,peak_loc_GW:length(signal_separated(1,:))),'r');
    title('Acceleration of rigth shank');
    xlabel('Samples');
    ylabel('Acc');

% -------------------------------------------------------------------------
% Interpolation of forceplate signal.
% -------------------------------------------------------------------------

% Resample the timesries on the time tange where the two time vectors
% overlap.
% ts_forceplate=timeseries(1:length(AP_forceplate)-1,...
%    AP_forceplate(1:length(AP_forceplate)-1));
% ts_GW=timeseries(1:length(signal_separated(1,:)),signal_separated(1,:));
% AP_FP_interpolated=synchronize(ts_forceplate,ts_GW,'Union');

% In this case, we are going to interpolate with the most duration signal.
% Set the vector with the points that we need to interpolate.
%div=peak_loc_FP/peak_loc_GW;

% Compute a linear interpolation between the begining and the peak
% location.
%AP_FP1=

% Campute a linear interoplation between the peak location and the end.
%div=(length(AP_forceplate)-peak_loc_FP)/(length(signal_separated(1,:))-peak_loc_GW);

%AP_FP2=

% Join both interpolations

%AP_FP_interpolated=[AP_FP1 AP_FP2];

% -------------------------------------------------------------------------
% Plot the first GW signal and forceplate signal
% -------------------------------------------------------------------------

% figure (3)
% subplot(3,1,1)
%     plot(AP_FP_interpolated);
%     title('Ant-Post COP excursions');
%     xlabel('Samples');
%     ylabel('AP-COP (mm)');
%     
% subplot(3,1,2)
%     plot(signal_separated(1,:),'r');
%     title('Pitch of rigth shank');
%     xlabel('Samples');
%     ylabel('Pitch (deg)');
% 
% subplot(3,1,3)
%     plot(AP_FP_interpolated)
%     hold on
%     plot(signal_separated(1,:),'r')
%     legend('AP excursions (forceplatform)','pitch rigth shank (GaitWatch)');
%     title('Synchronized Signals (GW and forceplatform)');
%     hold off
% 



% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END 
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\