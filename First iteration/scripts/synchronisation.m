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
% Last modification: 12/01/2015.
%
% -------------------------------------------------------------------------
%
% The present file synchronises the two data sets of GaitWatch and force
% plate, respectively. The file is structured as follows:
% 
% * 1) The GaitWatch library containing all the core functions is imported.
% 
% * 2) Load data from both force plate and GaitWatch. 
%
% * 3) Separate the ten cycles of the GaitWatch data. That is, find the
%      first peak of the ten cycles in the GaitWatch data and store the
%      following values to the first peak of the next cycle in a separate
%      vector.
% 
% * 4) Store the separate cycles in timeseries objects.
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

% Load the data from ES39_05 (forceplate) and GW_1605_150814_1154 (GW)
% We use pitch of rigth shank signal from GW and both feet AP from 
% foceplate to synchronize.


% force plate data in anterior-posterior direction. Change / against \ for 
% MATLAB on Windows!!!
load ('forcePlate_GW_data_ES39/both_AP.mat');                           
AP_FP=both_AP;                                  
load ('forcePlate_GW_data_ES39/timeFP');    
time_FP=timeFP;

% 
load ('forcePlate_GW_data_ES39/pitch_GKF_right_shank');
acc_GW=pitch_GKF_right_shank;               
load ('forcePlate_GW_data_ES39/time');
time_GW=time;

% -------------------------------------------------------------------------
% 3) Separate the ten cycles of the GaitWatch data. That is, find the
%    first peak of the ten cycles in the GaitWatch data and store the
%    following values to the first peak of the next cycle in a separate
%    vector. 
% -------------------------------------------------------------------------

% Find all peaks greater than 100 and store the first one in peak_loc_FP.
[peak_value_FP, peak_location_FP] = findpeaks(AP_FP, 'minpeakheight', 100);
peak_loc_FP=peak_location_FP(1);

% Select cycles manually by selecting points between two cycles.
indexes = gw.getDCindexes(acc_GW, 'SELECT POINTS BETWEEN CYCLES');
close all;

n_cycles = length(indexes)-1;

cycle_lengths = indexes(2:length(indexes))-indexes(1:length(indexes)-1);

GW_cycles = zeros(n_cycles, max(cycle_lengths));

for i = 1 : n_cycles
    
    % Find all peaks greater than 20 in cycle i.
    [peak_value_GW, peak_location_GW] = findpeaks(acc_GW, 'minpeakheight', 20);
    
    GW_cycles(i, : ) = acc_GW(peak_location_GW(1):indexes(i+1));
    
end

% -------------------------------------------------------------------------
% * 4) Store the separate cycles in timeseries objects.
% -------------------------------------------------------------------------

ts1_GW = timeseries(GW_cycles(1,:), 0:5:5*(max(size(GW_cycles))-1), 'name', 'GW_Cycle_1');
% ts2_GW = timeseries(GW_cycles(2,:), 0:5:5*(max(size(GW_cycles))-1));
% ts3_GW = timeseries(GW_cycles(3,:), 0:5:5*(max(size(GW_cycles))-1));
% ts4_GW = timeseries(GW_cycles(4,:), 0:5:5*(max(size(GW_cycles))-1));
% ts5_GW = timeseries(GW_cycles(5,:), 0:5:5*(max(size(GW_cycles))-1));
% ts6_GW = timeseries(GW_cycles(6,:), 0:5:5*(max(size(GW_cycles))-1));
% ts7_GW = timeseries(GW_cycles(7,:), 0:5:5*(max(size(GW_cycles))-1));
% ts8_GW = timeseries(GW_cycles(8,:), 0:5:5*(max(size(GW_cycles))-1));
% ts9_GW = timeseries(GW_cycles(9,:), 0:5:5*(max(size(GW_cycles))-1));
% ts10_GW = timeseries(GW_cycles(10,:), 0:5:5*(max(size(GW_cycles))-1));
% ts11_GW = timeseries(GW_cycles(11,:), 0:5:5*(max(size(GW_cycles))-1));
% ts12_GW = timeseries(GW_cycles(12,:), 0:5:5*(max(size(GW_cycles))-1));


    
   
    
     

   

% Create three timeseries objects to store the data collected at each intersection

count1 = timeseries(count(:,1), 1:24,'name', 'intersection1');
count2 = timeseries(count(:,2), 1:24,'name', 'intersection2');
count3 = timeseries(count(:,3), 1:24,'name', 'intersection3');

%open('count1')  % Display timeseries

%get(count1)     % Display properties

count1.DataInfo

count1.DataInfo.Units = 'cars';     % Change Units to cars

count1.DataInfo.Interpolation = tsdata.interpolation('zoh');    % Change interpolation to zero-order hold

count1.DataInfo

count1.TimeInfo.Units = 'hours';        % Modify the time units to be 'hours' for the three time series
count2.TimeInfo.Units = 'hours';
count3.TimeInfo.Units = 'hours';

% Add two events to the data that mark the times of the AM commute and PM commute.
% 
% Construct and add the first event to all time series. The first event occurs at 8 AM.

e1 = tsdata.event('AMCommute',8);
e1.Units = 'hours';            % Specify the units for time
count1 = addevent(count1,e1);  % Add the event to count1
count2 = addevent(count2,e1);  % Add the event to count2
count3 = addevent(count3,e1);  % Add the event to count3

% % Construct and add the second event to all time series. The second event occurs at 6 PM.
% 
e2 = tsdata.event('PMCommute',18);
e2.Units = 'hours';            % Specify the  units for time
count1 = addevent(count1,e2);  % Add the event to count1
count2 = addevent(count2,e2);  % Add the event to count2
count3 = addevent(count3,e2);  % Add the event to count3

% figure;
% % plot(count1);
% % hold on;
% % 
%  plot(count2, 'color', 'red');
% % 
% 
% plot(count3, 'color', 'green');

% Create a tscollection object nameed count_coll and use the constructor syntax

tsc = tscollection({count1 count2},'name', 'count_coll');

% Add the third timeseries object

tsc = addts(tsc, count3);

% Interpolate values at each half-hour mark.

tsc1 = resample(tsc,1:0.5:24);

% Plot the members of tsc1 with markers to see the results of interpolating.

plot(tsc1.intersection1,'-xb','Displayname','Intersection 1')
hold on
plot(tsc1.intersection2,'-.xm','Displayname','Intersection 2')
plot(tsc1.intersection3,':xr','Displayname','Intersection 3')
legend('show','Location','NorthWest')


% ------------- Synchronisation -----------------

% Create two timeseries, such that ts1.timeinfo.StartDate 
% is one day after ts2.timeinfo.StartDate:

ts1 = timeseries([1 2],[datestr(now); datestr(now+1)]);
ts2 = timeseries([1 2],[datestr(now-1); datestr(now)]);

% If you use this code, then ts1.timeinfo.StartDate
% is changed to match ts2.TimeInfo.StartDate 
% and ts1.Time changes to 1:

[ts1 ts2] = synchronize(ts1,ts2,'union');

% But if you use this code, then ts1.timeinfo.StartDate
% is unchanged and ts1.Time is still 0:

[ts1 ts2] = synchronize(ts1,ts2,'union','KeepOriginalTimes',true);



