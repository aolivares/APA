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
% Last modification: 14/01/2015.
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
% * 3) Separate the ten cycles of the GaitWatch data. That is, find the
%      first peak of the ten cycles in the GaitWatch data and store the
%      following values to the first peak of the next cycle in a separate
%      vector.
% 
% * 4) Store the separate cycles in a time series collection.
% 
% * 5) Resample the ten GaitWatch and the corresponding ten force plate
%      timeseries with a common time vector. 
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

fs_FP = 120;
fs_GW = 200;


% -------------------------------------------------------------------------
% 3) Select the cycles in the GaitWatch data manually and find peaks,
%    that is, points when patient touches the force plate. 
% -------------------------------------------------------------------------

% Select the cycles in the GaitWatch data manually by selecting points
% before and after each cycle.
indexes = gw.getDCindexes(a_Z_right_shank_1_C, ... 
          'Acceleration Z-axis right shank - Press ALT-key and select cycles');
      
n_cycles = length(indexes)-1;
peak_ind = zeros(1, n_cycles);

for i = 1 : n_cycles
    
    % Find all peaks greater than 1.4 in cycle i and store them in peak_ind.
    [peak_values, peak_locations] = findpeaks(a_Z_right_shank_1_C(indexes(i):indexes(i+1)), ...
                                    'minpeakheight', 1.4);
    
    peak_ind(i) = peak_locations(1)+indexes(i);
    
end

% Plot detected peaks
plot(time, a_Z_right_shank_1_C);
hold on;
plot(time(peak_ind), a_Z_right_shank_1_C(peak_ind), 'r.');


% -------------------------------------------------------------------------
% * 4) Store GaitWatch Data and force plate data in time series collections
% -------------------------------------------------------------------------

% Create time series of acceleration trunk.
a_trunk = timeseries([a_X_center_trunk_3_C; a_Y_center_trunk_3_C; ...
                      a_Z_center_trunk_3_C], time, ...
                      'name', 'Acceleration trunk');
a_trunk.TimeInfo.Units = 'seconds';
a_trunk.DataInfo.Units = 'g';

% Add events (touch of force plate) to time series acceleration trunk.
for i = 1:length(peak_ind)
    
    event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
            time(peak_ind(i)));
        
    event.Units = 'seconds';
    a_trunk = addevent(a_trunk, event);

end


% Create time series of acceleration thigh.
a_thigh = timeseries([a_X_left_thigh_1_C';  a_Z_left_thigh_1_C'; ...
                      a_X_right_thigh_1_C'; a_Z_right_thigh_1_C'], time, ...
                      'name', 'Acceleration thigh');
a_thigh.TimeInfo.Units = 'seconds';
a_thigh.DataInfo.Units = 'g';

% Add events (touch of force plate) to time series acceleration thigh.
for i = 1:length(peak_ind)
    
    event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
            time(peak_ind(i)));
        
    event.Units = 'seconds';
    a_thigh = addevent(a_thigh, event);

end


% Create time series of acceleration shank.
a_shank = timeseries([a_X_left_shank_1_C';  a_Z_left_shank_1_C'; ...
                      a_X_right_shank_1_C'; a_Z_right_shank_1_C'], time, ...
                      'name', 'Acceleration shank');
a_shank.TimeInfo.Units = 'seconds';
a_shank.DataInfo.Units = 'g';

% Add events (touch of force plate) to time series acceleration shank.
for i = 1:length(peak_ind)
    
    event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
            time(peak_ind(i)));
        
    event.Units = 'seconds';
    a_shank = addevent(a_shank, event);

end

% Store all accelerations in time series collection.
a_tsc = tscollection({a_trunk, a_thigh, a_shank});

hold off;
figure();
plot(a_trunk);


% Create empty time-series object and name it.
data_FP_concat_rs = timeseries();

for i = 1:n_cycles

% Extract force plate time vector from cycle i of the force plate data
% and scale it from milliseconds to seconds.
time_FP = (FP_data_complete{i, 1}(FP_data_complete{i, 2}:FP_data_complete{i, 3}, 1)')/1000;

% Calculate bias and correct time axis of force plate data.
time_bias = time(peak_ind(i))-time_FP(1);

time_FP_corr = time_FP + time_bias;

% Create time series of time-corrected first force plate cycle.
data_FP = timeseries(FP_data_complete{i, 1}(FP_data_complete{i, 2}:FP_data_complete{i, 3}, 2:5)', ...
                     time_FP_corr , 'name', strcat(num2str(i), '. force plate cycle'));
data_FP.TimeInfo.Units = 'seconds';
data_FP.DataInfo.Units = 'N';

% Add event (point in time when patient touches the force plate)
event = tsdata.event(strcat(num2str(i), '. touch of force plate'), time(peak_ind(i)));
event.Units = 'seconds';
data_FP = addevent(data_FP, event);

if i == 1               % Only for verification
    figure();
    plot(data_FP)        
end

% Resample force plate cycle with GaitWatch time axis.
data_FP_rs = resample(data_FP, time_FP_corr:1/fs_GW:time_FP_corr(length(time_FP_corr)));

% Concatenate cycle i with the previous force plate cycles.
data_FP_concat_rs = append(data_FP_concat_rs, data_FP_rs);

if i == 1               % Only for verification
    figure();
    plot(data_FP_concat_rs)        
end
   

end

set(data_FP_concat_rs, 'name', 'Force sensor data synchronised');

figure();
plot(data_FP_concat_rs);

