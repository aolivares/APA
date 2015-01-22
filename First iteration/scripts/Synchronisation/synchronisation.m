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
% Version:  2.0
%
% Last modification: 22/01/2015.
%
% -------------------------------------------------------------------------
%
% The present file synchronises the two data sets of GaitWatch and force
% plate. The file is structured as follows:
% 
% * 1) Import the GaitWatch library containing all the core functions.
% 
% * 2) Load data from both force plate and GaitWatch. 
%
% * 3) Find the first peak in each cycle, that is, the point in which the
%      patient walks on the force plate. 
% 
% * 4) Store GaitWatch signals and force plate signals in time series
%      collections.
% 
% * 5) Synchronise and store seperate force plate cycles in a time series
%      object.  
% 
% -------------------------------------------------------------------------
% 0) Clear workspace and close all figures.
% -------------------------------------------------------------------------

clear all; close all; clc;

% -------------------------------------------------------------------------
% 1) Import GaitWatch functions library.
% -------------------------------------------------------------------------
% From now on, all the functions have to be called using 'gw.functionName'.

gw = gwLibrary;

% -------------------------------------------------------------------------
% 2) Load data from both force plate and GaitWatch.
% -------------------------------------------------------------------------

% Select force plate data file with a dialog box (only .mat files).
[filename_FP, filepath_FP] = uigetfile('*.mat', ...
    'Select Force Plate data file (.mat)', '../../data/ForcePlate');

% Select GaitWatch data file with a dialog box (only .mat files).
[filename_GW, filepath_GW] = uigetfile('*.mat', ...
    'Select GaitWatch data file (.mat)', '../../data/GaitWatch');

% Load force plate data into workspace.
load(fullfile(filepath_FP, filename_FP));

% Load GaitWatch data into workspace.
load(fullfile(filepath_GW, filename_GW));

% Set sampling frequencies of the Force Plate and the Gaitwatch.
fs_FP = 120;
fs_GW = 200;

% -------------------------------------------------------------------------
% 3) Find the first peak in each cycle, that is, the point in which the
%    patient walks on the force plate. 
% -------------------------------------------------------------------------

% Set threshold for peak detection in acceleration signal and number of
% samples between two cycles.
threshold = 1.3;
gap = 1200;

% Find all peaks greater than threshold.
[peak_values, peak_locations] = findpeaks(a_Z_right_shank_1_C, ...
                                'minpeakheight', threshold);
                         
% Compute distance between two peaks.                             
peak_distance = diff(peak_locations);

% Create logical vector to select last peak of each cycle. That is the one
% before a gap of at least 1200 samples, respectively. 
select = peak_distance > gap;

% Shift logical vector by one to the right to select the first instead of
% the last peak of a cycle.
first_peak = [0, select(1:length(select) - 1)'];

% Add the very first detected peak.
peaks = [peak_locations(1)', peak_locations(logical(first_peak))'];

% Plot for verification.
plot(time, a_Z_right_shank_1_C);
hold on;
plot(time(peaks), a_Z_right_shank_1_C(peaks), 'r.');
hold off;

% -------------------------------------------------------------------------
% 4) Store GaitWatch signals and force plate signals in time series
%    collections.
% -------------------------------------------------------------------------

% Create time series of acceleration trunk.
a_trunk = createTimeseries([a_X_center_trunk_3_C; a_Y_center_trunk_3_C; ...
                            a_Z_center_trunk_3_C], time, ...
                            'Acceleration trunk', 'seconds', 'g', ...
                            time(peaks), '. touch of force plate');
                        
% Create time series of acceleration thigh.
a_thigh = createTimeseries([a_X_left_thigh_1_C';  a_Z_left_thigh_1_C'; ...
                            a_X_right_thigh_1_C'; a_Z_right_thigh_1_C'], ...
                            time, 'Acceleration thigh', 'seconds', 'g', ...
                            time(peaks), '. touch of force plate');
                        
% Create time series of acceleration shank.
a_shank = createTimeseries([a_X_left_shank_1_C';  a_Z_left_shank_1_C'; ...
                            a_X_right_shank_1_C'; a_Z_right_shank_1_C'], ...
                            time, 'Acceleration shank', 'seconds', 'g', ...
                            time(peaks), '. touch of force plate');

% Store all accelerations in time series collection.
a_tsc = tscollection({a_trunk, a_thigh, a_shank});

hold off;
figure();
plot(a_trunk);

% -------------------------------------------------------------------------
% 5) Synchronise and store seperate force plate cycles in a time series
%    object.
% -------------------------------------------------------------------------

% Compute number of selected cycles.
n_cycles = length(peaks);

% Create empty time series object and name it.
data_FP_concat_rs = timeseries();

for i = 1:n_cycles

% Extract force plate time vector from cycle i of the force plate data
% and scale it from milliseconds to seconds.
time_FP = (FP_data_complete{i, 1}(FP_data_complete{i, 2}: ...
           FP_data_complete{i, 3}, 1)') / 1000;

% Calculate bias and correct time axis of force plate data.
time_bias = time(peaks(i)) - time_FP(1);

time_FP_corr = time_FP + time_bias;

% Create time series of time-corrected first force plate cycle.
data_FP = timeseries(FP_data_complete{i, 1}(FP_data_complete{i, 2}: ...
                 FP_data_complete{i, 3}, 2:5)', time_FP_corr , 'name', ...
                 strcat(num2str(i), '. force plate cycle'));
data_FP.TimeInfo.Units = 'seconds';
data_FP.DataInfo.Units = 'N';

% Add event (point in time when patient touches the force plate)
event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
                     time(peaks(i)));
event.Units = 'seconds';
data_FP = addevent(data_FP, event);

if i == 1               % Only for verification
    figure();
    plot(data_FP)        
end

% Resample force plate cycle with GaitWatch time axis.
data_FP_rs = resample(data_FP, time_FP_corr : 1 / fs_GW : ...
                      time_FP_corr(length(time_FP_corr)));

% Concatenate cycle i with the previous force plate cycles.
data_FP_concat_rs = append(data_FP_concat_rs, data_FP_rs);

if i == 1               % Only for verification
    figure();
    plot(data_FP_concat_rs)        
end
   
end

% Set name of the time series.
set(data_FP_concat_rs, 'name', 'Force sensor data synchronised');

figure();
plot(data_FP_concat_rs);

