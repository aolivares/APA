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
    'Select Force Plate data file (.mat)', '../data/ForcePlate');

% Select GaitWatch data file with a dialog box (only .mat files).
[filename_GW, filepath_GW] = uigetfile('*.mat', ...
    'Select GaitWatch data file (.mat)', '../data/GaitWatch');

% Load force plate data into workspace.
load(fullfile(filepath_FP, filename_FP));

% Load GaitWatch data into workspace.
load(fullfile(filepath_GW, filename_GW));

%% 
% -------------------------------------------------------------------------
% 3) Select the cycles in the GaitWatch data manually. 
% -------------------------------------------------------------------------

% Localise seperate cycles by manually selecting points between two cycles.
indexes = gw.getDCindexes(a_Z_right_shank_1, 'Acceleration Z-axis right shank');

%% 

n_cycles = length(indexes)-1;
peak_ind = zeros(1, n_cycles);

for i = 1 : n_cycles
    
    % Find all peaks greater than 1.4 in cycle i and store them in peak_ind.
    [peak_values, peak_locations] = findpeaks(a_Z_right_shank_1(indexes(i):indexes(i+1)), 'minpeakheight', 1.4);
    
    peak_ind(i) = peak_locations(1);
   
end

%% 

% -------------------------------------------------------------------------
% * 4) Store the separate cycles in time series collection.
% -------------------------------------------------------------------------


acc_trunk = timeseries([], time);

              



%%


%GW_cycles = cell(n_cycles, 1);
% Store data from the first peak to the beginning of the next cycle in
    % cell array.
    %GW_cycles(i, 1) = acc_GW(peak_locations_GW(1):indexes(i+1));
    
%     temp_cycle = acc_GW(peak_locations_GW(1):indexes(i+1));
%     
%     GW_cycles(i, : ) = [temp_cycle, zeros(1, (max(cycle_lengths)-length(temp_cycle)))];


% -------------------------------------------------------------------------
% * 4) Store the separate cycles in time series collection.
% -------------------------------------------------------------------------

GW_cycles_ts = cell(n_cycles, 1);

for i = 1 : n_cycles
    
    % Store the separate cycles in time series objects arranged in a cell
    % array
    GW_cycles_ts(i, 1) = timeseries(GW_cycles(i,1), 0:1/fs_GW:(1/fs_GW)*length(GW_cycles(i, 1)), 'name', strcat('GW_Cycle_', num2str(i)));
    ts1_GW.TimeInfo.Units = 'milliseconds';

end

% Store the separate cycles in time series collection.
GW_cycles_tsc = tscollection(GW_cycles_ts, 'name', strcat('GW_cycles_time_series_collection of all', num2str(n_cycles),'cycles'));


%% 

% -------------------------------------------------------------------------
% 5) Resample the ten GaitWatch and the corresponding ten force plate
%    timeseries with a common time vector. 
% -------------------------------------------------------------------------

%size(GW_cycles_tsc)

%common_time = 0:5:min(length());

GW_cycles_tsc_sync = resample(GW_cycles_tsc, common_time);

That is, find the
%    first peak of the ten cycles, respectively, in the GaitWatch data
%    and store the following values to the first peak of the next cycle
%    in separate vectors in a cell array.  


