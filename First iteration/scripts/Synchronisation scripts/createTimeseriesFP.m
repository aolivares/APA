function [ ts_cell ] = createTimeseriesFP(data_FP, time_FP, ...
                                          sync_peak_times, fs_GW, name, ...
                                          time_unit, data_unit)
%CREATETIMESERIES creates a time series with the desired time and data
%units, and adds events.
%
% - Input:
%    |_ 'data_FP': Nx1 cell array consisting of N FP cycles.
%    |_ 'time_FP': Row vector with the length M that contains the time steps.
%    |_ 'sync_peaks_times': Row vector that contains the times of the
%    synchronisation peaks.
%    |_ 'fs_GW': Sampling frequency of the GaitWatch.
%    |_ 'name': String that serves as the name of the time series.
%    |_ 'time_unit': String containig the unit of the time vector.
%    |_ 'data_unit': String containig the unit of the data matrix.
%    
% - Output:
%    |_ 'ts': Cell array of synchronised and resampled time series objects
%    of the separate cycles containing the data and time vector as well as
%    units and events.  
%
% -------------------------------------------------------------------------
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Veronica Torres (2): 
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
% * Version:  2.1
%
% * Last modification: 23/01/2015
% -------------------------------------------------------------------------

% Compute number of cycles.
n_cycles = max(size(time_FP));

% Create empty cell array.
ts_cell = cell(1, n_cycles);

for i = 1:n_cycles

% Calculate bias to correct time vector of force plate data.
time_bias = sync_peak_times(i) - time_FP{i, 1}(1);

time_FP_corr = time_FP{i, 1} + time_bias;

% Create time series object
ts = timeseries(data_FP{i, 1}, time_FP_corr, 'name', strcat(num2str(i), ...
                ['. cycle of ', name]));
ts.TimeInfo.Units = time_unit;
ts.DataInfo.Units = data_unit;

% Add events to time series
event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
                     sync_peak_times(i));

event.Units = time_unit;
ts = addevent(ts, event);

% Resample time series with GaitWatch time vector.
ts_rs = resample(ts, time_FP_corr : 1/fs_GW : ...
                 time_FP_corr(length(time_FP_corr)));

% Add time series to cell array.             
ts_cell{1, i} = ts_rs;

end

