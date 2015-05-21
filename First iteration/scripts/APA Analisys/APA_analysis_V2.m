% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------------- APA Analysis --------------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accele,rometric Data.
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
% * Last modification: 06/05/2015

% INFORMATION: This file contains the routine to detect when the second
% step happens, determine the APAs of the FP and GW signals in this case 
% and the correlation between them. 

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

% Set flags which control the visibility of the figures.
showPlotsCheck = 'yes';
showPlotsAPA = 'yes';
showPlotsCorr = 'no';

% -------------------------------------------------------------------------
% 1) Select the .mat file and extrat the data form timeseries.
% -------------------------------------------------------------------------

% Selection and load of the synchronised file.
[filename, filepath] = uigetfile('*.mat', ...
    'Select the data file from the patient (.mat)', '../../data/Synchronised/');

load(fullfile(filepath,filename));

% Extract data from timeseries for plot.

% Sum of force of all force plate cells.
force_sum_complete_ts = append(force_sum_ts{1, :});
force_sum_data = force_sum_complete_ts.data;

% Force of right-left and front-back feet (four signals).
force_sensors_complete_ts = append(force_sensors_ts{1, :});
fs_data = force_sensors_complete_ts.data;

% Antereo-Posterior center of pressure.
AP_COP_complete_ts = append(AP_COP_ts{1, :});
AP_COP_data = AP_COP_complete_ts.data;

% Medio-Lateral center of pressure.
ML_COP_complete_ts = append(ML_COP_ts{1, :});
ML_COP_data = ML_COP_complete_ts.data;

% Acceleration of the shanks (X-Z axes and left-right direction).
a_shanks_complete_ts = append(a_shanks{1, :});
a_shanks_data = a_shanks_complete_ts.data;

% Acceleration of the trunk (X, Y, Z axes).
a_trunk_complete_ts = append(a_trunk{1, :});
a_trunk_data = a_trunk_complete_ts.data;

% Angular Velocity of trunk (x, Y, Z axes).
g_trunk_complete_ts = append(g_trunk{1, :});
g_trunk_data = g_trunk_complete_ts.data;

% -------------------------------------------------------------------------
% 2) Determine when the second step occurred.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 2.1) Determine the edges (interval) for the second step.
% -------------------------------------------------------------------------

% Determine when the second step occurred. We extract the part of the signal
% when the patient carried out the step. It happens between the beginning of
% activity period of each cycle and the end of FP data.
% We use the edges of the activity detection in the left shank because
% usually used to be clearer. But the difference between edges of the right
% and left shank is very small.

% We consider the second step starts in the beginning of the second period
% of activity of each cycle.
init_second_step = edges_left(3:4:length(edges_left));

% Determine the time point where there isn't force signal (considerating 
% both feet) for each cycle.
force_sum_data = reshape(force_sum_data(3,1,:), 1, ...
                 length(force_sum_data(3,1,:)));

edges_FP_signal = find(diff(force_sum_data > 100)~=0);
final_second_step_edge = edges_FP_signal(2:2:...
                        length(edges_FP_signal));
final_second_step = force_sum_complete_ts.time(final_second_step_edge);

%--------------------------------------------------------------------------
% 3) Detect APA in trunk signal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 3.1) Determine when the second step (when patient goes down from
% plateform) starts with left or right foot.
%--------------------------------------------------------------------------
% When the patient starts to step with left foot, last value of the cycle
% of ML COP is positive because the step finishes the pressure is located
% in the right foot.
% When the patient starts with right foot, the last value of the ML COP
% cycle will be negative.
AP_COP_data =  reshape(AP_COP_data(3, 1, :), ...
                [1, max(size(AP_COP_data))]);
            
ML_COP_data =  reshape(ML_COP_data(3, 1, :), ...
                [1, max(size(ML_COP_data))]);
            
cycle_start_right = find(ML_COP_data(final_second_step_edge) < 0);


%--------------------------------------------------------------------------
% 3.2) Determine the APA peaks in the Center of Pressure.
%--------------------------------------------------------------------------

% Firstly, we detect the APA peaks in Medio-Lateral direction because we
% only have one peak, so is easier to not confuse with the point of time
% when this happens. This point is exactly the same that in the AP
% direction. Therefore, we can use the result to calculate this after.

% If the patient starts to walk with the right foot, the APA peak is
% positive. However, if the patient starts to walk with the left foot,
% the APA peak is negative.

% Determine the initial and end point (time) of each interval where 
% we need to find the peaks, i.e, the second activity period of each cycle.
initcross = time_GW(init_second_step);
finalcross = final_second_step';

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point roundoff error.
    initcross_COP = find(abs(ML_COP_complete_ts.time - initcross(k))...
                    < 0.001);
    finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k))...
                    < 0.001);
                
     % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a positive peak.
        
         % Find all peaks in each interval.
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                ML_COP_data(...
                                initcross_COP:finalcross_COP));

        % Store the index of the longest positive peak.                                      
        peaks_APA_ML_COP(k) = find(ML_COP_data(...
                        initcross_COP:finalcross_COP)== max(...
                        pos_peak_values), 1) + initcross_COP - 1;
                    
    else % Look for a negative peak.
        
         % Find all peaks in each interval.    
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -ML_COP_data(...
                                initcross_COP:finalcross_COP));

        % Store the index of the longest negative peak.                                      
        peaks_APA_ML_COP(k) = find(ML_COP_data(...
                        initcross_COP:finalcross_COP)== -max(...
                        neg_peak_values), 1) + initcross_COP - 1;
                           
    end
end

% We calculate the value of the APA peaks.
value_APA_ML_COP = ML_COP_data(peaks_APA_ML_COP);
value_APA_AP_COP = AP_COP_data(peaks_APA_ML_COP);

%------------------------------- Plots-------------------------------------
if strcmpi(showPlotsCheck,'yes')
subplot(2,1,1)
plot(AP_COP_complete_ts.time, AP_COP_data, 'g');
hold on;
plot(AP_COP_complete_ts.time(peaks_APA_ML_COP), value_APA_AP_COP , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['AP COP with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('AP COP (mm)');

subplot(2,1,2)
plot(ML_COP_complete_ts.time, ML_COP_data, 'g');
hold on;
plot(ML_COP_complete_ts.time(peaks_APA_ML_COP), value_APA_ML_COP , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['ML COP with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('ML COP (mm)');
end

%--------------------------------------------------------------------------
% 3.3) Determine the next peak after the APA peaks detected in ML direction. 
% This point will indicate the interval to find the APA peaks in the 
% acceleration signals 
%--------------------------------------------------------------------------
initcross_COP = peaks_APA_ML_COP;

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point round off error.
    finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k))...
                    < 0.001);
                
     % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a negative peak.
        
         % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -ML_COP_data(...
                                initcross_COP(k):finalcross_COP));

        % Store the index of the first negative peak.
        neg_peak_locations = sort(neg_peak_locations);
        finalcross_acc(k) = neg_peak_locations(1) + initcross_COP(k) -1;

                    
    else % Look for a positive peak.
        
         % Find all peaks in each interval.
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                ML_COP_data(...
                                initcross_COP(k):finalcross_COP));

        % Store the index of the first negative peak.
        pos_peak_locations = sort(pos_peak_locations);
        finalcross_acc(k) = pos_peak_locations(1) + initcross_COP(k) -1;
                           
    end
end

                    
%--------------------------------------------------------------------------
% 3.4) Determine the APA peaks in the Acceleration Signals.
%--------------------------------------------------------------------------
% Calculate the indexes to be able to use them for the accelerations
% signals.
initcross_acc_time = ML_COP_complete_ts.time (initcross_COP);
finalcross_acc_time = ML_COP_complete_ts.time (finalcross_acc);

a_trunk_data_X = reshape(a_trunk_data(1, 1, :), ...
                [1, max(size(a_trunk_data))]);
a_trunk_data_Y = reshape(a_trunk_data(2, 1, :), ...
                [1, max(size(a_trunk_data))]);
            
for k = 1:length(initcross)
    
  % We obtain the interval where the APA peaks in the Acc signals appear.
  initcross_acc = find(abs( a_trunk_complete_ts.time- initcross_acc_time(k))...
                    < 0.001);
  finalcross_acc = find(abs( a_trunk_complete_ts.time- finalcross_acc_time(k))...
                    < 0.001);
                
  % We ontain the value of the peak in the AP direction. The position of
  % the peak is the same in the ML direction as well.
  % Find all peaks in each interval.    
    [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -a_trunk_data_X(...
                            initcross_acc:finalcross_acc));

   % Store the index of the longest negative peak.                                      
    peaks_APA_acc_X(k) = find(a_trunk_data_X(...
                    initcross_acc:finalcross_acc)== -max(...
                    neg_peak_values), 1) + initcross_acc - 1;
  
end

% We calculate the value of the APA peaks.
value_APA_acc_X = a_trunk_data_X(peaks_APA_acc_X);
value_APA_acc_Y = a_trunk_data_Y(peaks_APA_acc_X);

%------------------------------- Plots-------------------------------------
if strcmpi(showPlotsCheck,'yes')
figure ()
subplot(2,1,1)
plot(a_trunk_complete_ts.time, a_trunk_data_X, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_acc_X), value_APA_acc_X , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time and the APA peaks']);
xlabel('Time in s');
ylabel('Acceleration (g)');

subplot(2,1,2)
plot(a_trunk_complete_ts.time, a_trunk_data_Y, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_acc_X), value_APA_acc_Y , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Acceleration in (g)');
end

%--------------------------------------------------------------------------
% 3.5) Determine the APA peaks in the Gyroscope Signals.
%--------------------------------------------------------------------------
g_trunk_data_X = reshape(g_trunk_data(1, 1, :), ...
                [1, max(size(g_trunk_data))]);
g_trunk_data_Y = reshape(g_trunk_data(2, 1, :), ...
                [1, max(size(g_trunk_data))]);
            
            
            