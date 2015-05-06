function [ ts_cell ] = createTimeseriesGW(data_GW, time_GW, sync_peaks, ...
                                          last_peaks, add_samples, ...
                                          name, time_unit, data_unit)
%CREATETIMESERIES creates a cell array of time series of the seperate
%cycles of the GaitWatch signal swith the desired time and data units, and 
%adds the time of the synchronisation peaks as events.
%
% - Input:
%    |_ 'data_GW': NxM Data matrix consisting of N signals with M samples.
%    |_ 'time_GW': Row vector with the length M that contains the time steps.
%    |_ 'sync_peaks': Row vector that contains the indexes of the
%    synchronisation peaks.
%    |_ 'last_peaks': Row vector that contains the indexes of the
%    last peak of each cycle.
%    |_ 'add_samples': Positive integer that sets the number of samples
%    before sync_peaks and after last_peaks, respectively.
%    |_ 'name': String that serves as the name of the time series.
%    |_ 'time_unit': String containig the unit of the time vector.
%    |_ 'data_unit': String containig the unit of the data matrix.
%    
% - Output:
%    |_ 'ts': Cell array of time series objects of the separate cycles
%    containing the data and time vector as well as units and events.  
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
% * Last modification: 25/01/2015
% -------------------------------------------------------------------------

% Compute number of cycles.
n_cycles = length(sync_peaks);

% Create empty cell array.
ts_cell = cell(1, n_cycles);

for i = 1:n_cycles

% Create time series object 

  if( i == n_cycles)
    ts = timeseries(data_GW(:, sync_peaks(i)-add_samples:last_peaks(i)), ...
                time_GW(sync_peaks(i)-add_samples:last_peaks(i)), ...
                'name', strcat(num2str(i), ['. cycle of ', name]));
  else
    ts = timeseries(data_GW(:, sync_peaks(i)-add_samples:last_peaks(i)+add_samples), ...
                time_GW(sync_peaks(i)-add_samples:last_peaks(i)+add_samples), ...
                'name', strcat(num2str(i), ['. cycle of ', name]));  
  end
  
ts.TimeInfo.Units = time_unit;
ts.DataInfo.Units = data_unit;

% Add events to time series
event = tsdata.event(strcat(num2str(i), '. touch of force plate'), ...
                     time_GW(sync_peaks(i)));

event.Units = time_unit;
ts = addevent(ts, event);

% Add time series to cell array.
ts_cell{1, i} = ts;

end

