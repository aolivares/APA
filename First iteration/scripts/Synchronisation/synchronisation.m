% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ||||||||| Synchronisation of GaitWatch and force plate data |||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%
% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Veronica  Torres (2): 
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
% Version:  2.2
%
% Last modification: 25/01/2015.
%
% -------------------------------------------------------------------------
%
% The present file synchronises the two data sets of the GaitWatch and the
% force plate. The file is structured as follows:
% 
% * 1) Import the GaitWatch library containing all the core functions.
% 
% * 2) Load data from both the force plate and the GaitWatch. 
%
% * 3) Find the first peak of each cycle in the GaitWatch signal of the
%      acceleration of the shank, that is, the point in which the patient 
%      walks on the force plate. Detect if the patient steps with the left or
%      the right foot first and use this point for synchronisation. 
% 
% * 4) Store seperate GaitWatch cycles in a cell array of time series
%      objects and add the point in time when the patient walks on the
%      force plate as an event.
% 
% * 5) Synchronise and resample the seperate force plate cycles, then store 
%      them in a cell array of time series objects and add the point in
%      time when the patient walks on the force plate as an event. 
% 
% -------------------------------------------------------------------------
% 0) Clear workspace and close all figures.
% -------------------------------------------------------------------------

clear all; close all; clc;

% -------------------------------------------------------------------------
% 1) Import the GaitWatch library containing all the core functions.
% -------------------------------------------------------------------------
% From now on, all the functions have to be called using 'gw.functionName'.

gw = gwLibrary;

% -------------------------------------------------------------------------
% 2) Load data from both the force plate and the GaitWatch.
% -------------------------------------------------------------------------

% Select force plate data file with a dialog box (only .mat files).
[filename_FP, filepath_FP] = uigetfile('*.mat', ...
'Select Force Plate data file (.mat)', ...
'/Users/Rob/Documents/APA/First iteration/data/ForcePlate/Preprocessed');

% Select GaitWatch data file with a dialog box (only .mat files).
[filename_GW, filepath_GW] = uigetfile('*.mat', ...
'Select GaitWatch data file (.mat)', ...
'/Users/Rob/Documents/APA/First iteration/data/GaitWatch/Calibrated');

% Load force plate data into workspace.
load(fullfile(filepath_FP, filename_FP));

% Load GaitWatch data into workspace.
load(fullfile(filepath_GW, filename_GW));

% Set sampling frequencies of the Force Plate and the GaitWatch.
fs_FP = 120;
fs_GW = 200;

% Signal mapping for testing.
time_GW = time;     

% -------------------------------------------------------------------------
% 3) Find the first peak of each cycle in the GaitWatch signal of the
%    acceleration of the shank, that is, the point in which the patient 
%    walks on the force plate. Detect if the patient steps with the left or
%    the right foot first and use this point for synchronisation.
% -------------------------------------------------------------------------

% Set threshold for peak detection in acceleration signal and minimum
% number of samples between two cycles.
threshold = 1.3;
gap = 1300;

% Find all peaks greater than threshold.
[peak_values_l, peak_locations_l] = findpeaks(a_Z_left_shank_1_C, ...
                                              'minpeakheight', threshold);
[peak_values_r, peak_locations_r] = findpeaks(a_Z_right_shank_1_C, ...
                                              'minpeakheight', threshold);

% Plot for verification.
close all;
plot(time, a_Z_right_shank_1_C);
hold on;
plot(time(peak_locations_r), a_Z_right_shank_1_C(peak_locations_r), 'r.')
figure();
plot(time, a_Z_left_shank_1_C, 'g');
hold on;
plot(time(peak_locations_l), a_Z_left_shank_1_C(peak_locations_l), 'k.');

%%

% Compute distance between two peaks.                             
peak_distance_l = diff(peak_locations_l);
peak_distance_r = diff(peak_locations_r);

% Create logical vector to select last peak of each cycle. That is the one
% before a gap of at least the number of samples stored in gap,
% respectively. 
select_last_l = peak_distance_l > gap;
select_last_r = peak_distance_r > gap;

last_peaks_l = peak_locations_l(select_last_l);
last_peaks_r = peak_locations_r(select_last_r);

close all;
plot(time, a_Z_right_shank_1_C);
hold on;
plot(time(last_peaks_r), a_Z_right_shank_1_C(last_peaks_r), 'r.');
figure()
plot(time, a_Z_left_shank_1_C, 'g');
hold on;
plot(time(last_peaks_l), a_Z_left_shank_1_C(last_peaks_l), 'k.');

%%

% Shift logical vector by one to the right to select the first instead of
% the last peak of a cycle.
select_first_l = [0, select_last_l(1:length(select_last_l) - 1)'];
select_first_r = [0, select_last_r(1:length(select_last_r) - 1)'];

% Add the very first detected peak.
sync_peaks_l = [peak_locations_l(1)', peak_locations_l(logical(select_first_l))'];
sync_peaks_r = [peak_locations_r(1)', peak_locations_r(logical(select_first_r))'];

% Plot for verification.
close all;
plot(time, a_Z_right_shank_1_C);
hold on;
plot(time(sync_peaks_r), a_Z_right_shank_1_C(sync_peaks_r), 'r.');
figure();
plot(time, a_Z_left_shank_1_C, 'g');
hold on;
plot(time(sync_peaks_l), a_Z_left_shank_1_C(sync_peaks_l), 'k.');

%%

% Evaluate if the patient steps with the left or right limb first for each
% cycle and store the first peak in sync_peaks, respectively, then sort it.
sync_peaks = [sync_peaks_r(sync_peaks_l > sync_peaks_r), ...
              sync_peaks_l(sync_peaks_r > sync_peaks_l)];
sync_peaks = sort(sync_peaks);

last_peaks = [last_peaks_r(last_peaks_l < last_peaks_r)', ...
              last_peaks_l(last_peaks_r < last_peaks_l)'];
last_peaks = sort(last_peaks);

% Plot for verification.
close all;
hold on;
plot(time, a_Z_left_shank_1_C);
plot(time(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)), 'r.');
plot(time, a_Z_right_shank_1_C, 'g');
plot(time(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)), 'k.');

%%

% -------------------------------------------------------------------------
% 4) Store seperate GaitWatch cycles in a cell array of time series
%    objects and add the point in time when the patient walks on the
%    force plate as an event.
% -------------------------------------------------------------------------

% Set the additional time in seconds before (sync_peak) and after
% (last_peak) each cycle that will be stored in the time series.
add_time = 2;
add_samples = floor(add_time * fs_GW);

% Create cell array containing the time series of the seperate cycles of 
% acceleration trunk.
a_trunk = createTimeseriesGW([a_X_center_trunk_3_C; a_Y_center_trunk_3_C; ...
                              a_Z_center_trunk_3_C], time, sync_peaks, ...
                              last_peaks, add_samples, ...
                              'Acceleration trunk', 'seconds', 'g');
                        
% Create cell array containing the time series of the seperate cycles of 
% acceleration thigh.
a_thigh = createTimeseriesGW([a_X_left_thigh_1_C';  a_Z_left_thigh_1_C'; ...
                              a_X_right_thigh_1_C'; a_Z_right_thigh_1_C'], ...
                              time, sync_peaks, last_peaks, add_samples, ...
                              'Acceleration thigh', 'seconds', 'g');
                        
% Create cell array containing the time series of the seperate cycles of 
% acceleration shank.
a_shank = createTimeseriesGW([a_X_left_shank_1_C';  a_Z_left_shank_1_C'; ...
                              a_X_right_shank_1_C'; a_Z_right_shank_1_C'], ...
                              time, sync_peaks, last_peaks, add_samples, ...
                              'Acceleration shank', 'seconds', 'g');

% Plot for verification.
close all;
plot(a_trunk{1, 1});

%%

% -------------------------------------------------------------------------
% 5) Synchronise and resample the seperate force plate cycles, then store 
%    them in a cell array of time series objects and add the point in
%    time when the patient walks on the force plate as an event.
% -------------------------------------------------------------------------

% Rearrange data for testing until data extraction script is finished.
% Compute number of cycles.
n_cycles = length(sync_peaks);

FP_data = cell(n_cycles, 1);
FP_time = cell(n_cycles, 1);

for i=1:n_cycles
    
    FP_data(i, 1) = {FP_data_complete{i, 1}(FP_data_complete{i, 2}: ...
                     FP_data_complete{i, 3}, 2:5)};
             
    FP_time(i, 1) = {FP_data_complete{i, 1}(FP_data_complete{i, 2}: ...
                     FP_data_complete{i, 3}, 1) / 1000};
    
end

%Plot for verification.
close all;
plot(FP_time{1, 1}, FP_data{1, 1});

%%

% Compute the times where the syncronisation peaks appear from the
% sync_peaks vector containing the indexes of the peaks and the GW time
% vector.
sync_peak_times = time_GW(sync_peaks);

% Create cell array containing the time series of the seperate cycles of
% the four force sensors.
force_sensors = createTimeseriesFP(FP_data, FP_time, sync_peak_times, ...
                                   fs_GW, 'Force Sensors', 'seconds', 'N');

% Plot for verification.
close all;
plot(force_sensors{1, 1});
figure();
plot(a_trunk{1, 1});

%%

% -------------------------------------------------------------------------
% 6) Store time series objects as .mat file
% -------------------------------------------------------------------------

save('/Users/Rob/Documents/APA/First iteration/data/GaitWatch/Synchronised/GW_sync_ES_39.mat', ...
     'a_trunk','a_shank', 'a_thigh', ...
     'force_sensors');