function [ ts ] = createTimeseries(data, time, name, time_unit, ...
                                   data_unit, event_time, event_name)
%CREATETIMESERIES creates a time series with the desired time and data
%units, and adds events.
%
% - Input:
%    |_ 'data': NxM Data matrix consisting of N signals with M samples.
%    |_ 'time': The time vector is a row vector with the length M.
% 
%    |_ 'time_unit': String containig the unit of the time vector.
%    |_ 'data_unit': String containig the unit of the data matrix.
%    |_ 'event_time': Row vector containing the event times.
%    |_ 'event_name': String containing the event name.
%    
% - Output:
%    |_ 'ts': Time series object containing the data and time vector as
%             well as units and events.  
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
% * Last modification: 19/01/2015
% -------------------------------------------------------------------------

% Create time series object
ts = timeseries(data, time, 'name', name);
ts.TimeInfo.Units = time_unit;
ts.DataInfo.Units = data_unit;

% Add events to time series
for i = 1:length(event_time)
    
    event = tsdata.event(strcat(num2str(i), event_name), event_time(i));
        
    event.Units = time_unit;
    ts = addevent(ts, event);

end


end

