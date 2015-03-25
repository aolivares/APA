% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------------- APA Analysis --------------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
% * Last modification: 20/03/2015

% INFORMATION:

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

[filename, filepath] = uigetfile('*.mat', ...
    'Select the data file from the patient (.mat)', '../../data/Synchronised/');

load(fullfile(filepath,filename));

% -------------------------------------------------------------------------
% 1) Determine APA characteristics in trunk signals, COP and the
% correlation between them.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 1.1) Determine when the second step ocurred
% -------------------------------------------------------------------------

% Determine when the second step ocurred. We extract the part of the signal
% when the patient carried out the step. It happens between the beginning of
% activity period of each cycle and the end of FP data.
% We use the edges of the activity detection in the right shank because
% usually used to be clearer. But the difference between edges of the right
% and left shank is very small.

% We consider the second step starts in the beginning of the second period
% of activity of each cycle.
init_second_step = edges_right(3:4:length(edges_right));

% Extract data from timeseries for plot.
force_sum_complete_ts = append(force_sum_ts{1, :});
force_sum_data = force_sum_complete_ts.data;

force_sensors_complete_ts = append(force_sensors_ts{1, :});
fs_data = force_sensors_complete_ts.data;

AP_COP_complete_ts = append(AP_COP_ts{1, :});
AP_COP_data = AP_COP_complete_ts.data;

ML_COP_complete_ts = append(ML_COP_ts{1, :});
ML_COP_data = ML_COP_complete_ts.data;

a_shanks_complete_ts = append(a_shanks{1, :});
a_shanks_data = a_shanks_complete_ts.data;

a_trunk_complete_ts = append(a_trunk{1, :});
a_trunk_data = a_trunk_complete_ts.data;

% Determine the time point where there isn't force signal (considerating 
% both feet) for each cycle.
force_sum_data = reshape(force_sum_data(3,1,:), 1, ...
                 length(force_sum_data(3,1,:)));

final_second_step_edge = find(diff(force_sum_data > 100)~=0);
final_second_step_edge = final_second_step_edge(2:2:...
                        length(final_second_step_edge));
final_second_step = force_sum_complete_ts.time(final_second_step_edge);


% -------------------------------------------------------------------------
% 1.2) Plots
% -------------------------------------------------------------------------

% --------------- Acceleration in shanks and force-------------------------

figure();
subplot(3, 1, 1);
plot(a_shanks_complete_ts.time, reshape(a_shanks_data(2, 1, :), ...
     [1, max(size(a_shanks_data))]));

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the left shank with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');

subplot(3, 1, 2)
plot(a_shanks_complete_ts.time, reshape(a_shanks_data(4, 1, :), ...
     [1, max(size(a_shanks_data))]), 'g');

 
% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the right shank with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');

subplot(3, 1, 3)
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('The force in FP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('Force in N');

% --------------------Acceleration in trunk and COP------------------------

figure()
subplot(4, 1, 1);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(1, 1, :), ...
     [1, max(size(a_trunk_data))]))

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');


subplot(4, 1, 2)
plot(AP_COP_complete_ts.time, reshape(AP_COP_data(3, 1, :), ...
     [1, max(size(AP_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('AP COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');

subplot(4, 1, 3);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(2, 1, :), ...
     [1, max(size(a_trunk_data))]));

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');

subplot(4, 1, 4)
plot(ML_COP_complete_ts.time, reshape(ML_COP_data(3, 1, :), ...
     [1, max(size(ML_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('ML COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');

% --------------------Interesting cycles (Zoom)----------------------------
% ------------------------> 6º cycle
figure();
subplot(5, 1, 1);
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('The force in FP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('Force in N');
axis([145, 170, 0, 1500]);

subplot(5, 1, 2);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(1, 1, :), ...
     [1, max(size(a_trunk_data))]))

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([145, 170, 0, 1]);


subplot(5, 1, 3)
plot(AP_COP_complete_ts.time, reshape(AP_COP_data(3, 1, :), ...
     [1, max(size(AP_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('AP COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');
axis([145, 170, 0, 400]);

subplot(5, 1, 4);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(2, 1, :), ...
     [1, max(size(a_trunk_data))]));

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([145, 170, -1, 1]);

subplot(5, 1, 5)
plot(ML_COP_complete_ts.time, reshape(ML_COP_data(3, 1, :), ...
     [1, max(size(ML_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('ML COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');
axis([145, 170, -200, 150]);

% ------------------------> 8º cycle
figure();
subplot(5, 1, 1);
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('The force in FP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('Force in N');
axis([190, 210, 0, 1500]);

subplot(5, 1, 2);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(1, 1, :), ...
     [1, max(size(a_trunk_data))]))

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([190, 210, 0, 1]);


subplot(5, 1, 3)
plot(AP_COP_complete_ts.time, reshape(AP_COP_data(3, 1, :), ...
     [1, max(size(AP_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('AP COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');
axis([190, 210, 0, 400]);

subplot(5, 1, 4);
plot(a_trunk_complete_ts.time, reshape(a_trunk_data(2, 1, :), ...
     [1, max(size(a_trunk_data))]));

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([190, 210, -1, 1]);

subplot(5, 1, 5)
plot(ML_COP_complete_ts.time, reshape(ML_COP_data(3, 1, :), ...
     [1, max(size(ML_COP_data))]));

% Vertical line at init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('ML COP with lines marker when the patient steps with the second time');
xlabel('Time in s');
ylabel('COP in mm');
axis([190, 210, -200, 150]);

% -------------------------------------------------------------------------
% 1.3) Differences when the patient starts with left or right foot in the 
% first step to see the differences.
% -------------------------------------------------------------------------

% Determine the cycles when patient start with left or right foot.

% Beginning the step with the left foot.
cycle_start_left = find(sync_peaks_l < sync_peaks_r);
a_trunk_left = a_trunk(cycle_start_left);

% Inicialization of the variables
a_trunk_left_X = zeros(length(cycle_start_left),5000); 
a_trunk_left_Y = zeros(length(cycle_start_left),5000);

% We group the cycles start with the left foot.
for i = 1:length(cycle_start_left)
    
   a_trunk_left_data = a_trunk_left{i}.data;
   
   a_trunk_left_data_X = reshape(a_trunk_left_data(1, 1, :), ...
       [1, max(size(a_trunk_left_data))]);
   a_trunk_left_X (i,1:length(a_trunk_left_data_X)) = a_trunk_left_data_X;
   
   a_trunk_left_data_Y = reshape(a_trunk_left_data(2, 1, :), ...
       [1, max(size(a_trunk_left_data))]);
   a_trunk_left_Y (i,1:length(a_trunk_left_data_Y)) = a_trunk_left_data_Y;  
    
end

figure()
plot(a_trunk_left_Y(2:3,:)');
title('Pattern of the Acceleration signal of the trunk in the Y-axe when step happend with the left foot');
axis([0, 2000, -0.5, 0.5]);


% Beginning with the right foot.
cycle_start_right = find(sync_peaks_r < sync_peaks_l);
a_trunk_right = a_trunk(cycle_start_right);


% Inicialization of the variables.
a_trunk_right_X = zeros(length(cycle_start_left),5000); 
a_trunk_right_Y = zeros(length(cycle_start_left),5000);

% We group the cycles start with the left foot.
for i = 1:length(cycle_start_right)
    
   a_trunk_right_data = a_trunk_right{i}.data;
   
   a_trunk_right_data_X = reshape(a_trunk_right_data(1, 1, :), ...
       [1, max(size(a_trunk_right_data))]);
   a_trunk_right_X (i,1:length(a_trunk_right_data_X)) = a_trunk_right_data_X;
   
   a_trunk_right_data_Y = reshape(a_trunk_right_data(2, 1, :), ...
       [1, max(size(a_trunk_right_data))]);
   a_trunk_right_Y (i,1:length(a_trunk_right_data_Y)) = a_trunk_right_data_Y;  
    
end

figure()
plot(a_trunk_right_Y(1:2,:)');
title('Pattern of the Acceleration signal of the trunk in the Y-axe when step happend with the right foot');
axis([0, 2000, -0.5, 0.5]);

figure ()
plot(a_trunk_right_X');
title('Pattern of the Acceleration signal of the trunk in the X-axe');
axis([0, 2000, 0, 1]);

% -------------------------------------------------------------------------
% 2) Detect APA in trunk signal.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 2.1) Application of Activity detectors in trunk signal
% -------------------------------------------------------------------------
% wag = wagLibrary;    
% 
% % Define input signal for rigth shank.
% axC =  reshape(a_trunk_data(1, 1, :), ...
%      [1, max(size(a_trunk_data))]);
% ayC =  reshape(a_trunk_data(2, 1, :), ...
%      [1, max(size(a_trunk_data))]);
% 
% input_signal = sqrt(axC .^ 2 + ayC .^ 2)';
% 
% % Computation of intensity markers.
% 
% % LTSD (window size, decision threshold and overlapping).
% lwin_ltsd = 100;       threshold_ltsd = 14;   shift_ltsd = 10;
% 
% % Get the decision signal of the LTSD algorithm and the marker.
% [V_ltsd, T_ltsd] = wag.ltsd(input_signal, lwin_ltsd, shift_ltsd, 512, ...
%     threshold_ltsd);
% [marker_ltsd, T_ltsd_expanded] = wag.compEstMark(V_ltsd, T_ltsd, ...
%     input_signal, lwin_ltsd, shift_ltsd);
%   
% figure
% plot(T_ltsd_expanded)
% hold on
% plot(threshold_ltsd * ones(1, length(T_ltsd_expanded)), 'r')
% legend('Detector output (LTSD)', 'Detection threshold')
% 
% figure
% plot(input_signal)
% hold on
% plot(marker_ltsd +1, 'r')
% legend('Input signal','LTSD decision')
% 
% % Determinate the initial and end point of each interval where we need to 
% % find the peaks, i.e, the first activity period of each cycle. 
% edges = find(diff(marker_ltsd)~=0);
% initcross = edges(3:4:length(edges));
% finalcross = edges(4:4:length(edges));

% -------------------------------------------------------------------------
% 2.2) Find the first 'negative' peak of each cycle in the signal 
%      of the x-axis of the acceleration of the trunk, that is, the point 
%      in time when the patient shifted backward. 
% -------------------------------------------------------------------------

% Determine the initial and end point (in time) of each interval where 
% we need to find the peaks, i.e, the second activity period of each cycle.
initcross = time_GW(init_second_step);
finalcross = final_second_step';

% Obtain the trunk signal of X-axis (Antereo-Posterior movement) and AP COP
% of both feet.
a_trunk_data_X = reshape(a_trunk_data(1, 1, :), ...
                [1, max(size(a_trunk_data))]);
AP_COP_data_shanks = reshape(AP_COP_data(3, 1, :), ...
                [1, max(size(AP_COP_data))]);

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point roundoff error.
    initcross_trunk = find(abs(a_trunk_complete_ts.time - initcross(k)) < 0.001);
    initcross_COP = find(abs(AP_COP_complete_ts.time - initcross(k)) < 0.001);

    finalcross_trunk = find(abs(a_trunk_complete_ts.time - finalcross(k)) < 0.001);
    finalcross_COP = find(abs(AP_COP_complete_ts.time - finalcross(k)) < 0.001);
    
    % Find all peaks in each interval.
    [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -a_trunk_data_X(...
                            initcross_trunk:finalcross_trunk));

    % Store the index of the first negative peak.                                      
    peaks_APA_trunk_X(k) = find(a_trunk_data_X(...
                    initcross_trunk:finalcross_trunk)== -max(...
                    neg_peak_values), 1, 'last') + initcross_trunk - 1;
                
     % Find all peaks in each interval.
    [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -AP_COP_data_shanks(...
                            initcross_COP:finalcross_COP));

    % Store the index of the first negative peak.                                      
    peaks_APA_AP_COP(k) = find(AP_COP_data_shanks(...
                    initcross_COP:finalcross_COP)== -max(...
                    neg_peak_values), 1, 'last') + initcross_COP - 1;
                 
end

% Calculate the value of the APA peaks in each cycle in the trunk 
% acceleration and the AP COP. 
value_APA_trunk_X = a_trunk_data_X(peaks_APA_trunk_X);
value_APA_AP_COP = AP_COP_data_shanks(peaks_APA_AP_COP);

% Calculate the linear correlation between the peaks detected with the 
% acceleration signal of the trunk and AP COP.
[corr_trunk_AP, prob_trunk_AP] = corr(value_APA_trunk_X', value_APA_AP_COP');

%--------------------------------------------------------------------------
% Plots
%--------------------------------------------------------------------------

%-------------Peaks of APA in trunk acceleration and AP COP----------------
figure ();
subplot(2,1,1)
plot(a_trunk_complete_ts.time, a_trunk_data_X, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_trunk_X), value_APA_trunk_X , 'r.');

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
ylabel('Acceleration (g)');


subplot(2,1,2)
plot(AP_COP_complete_ts.time, AP_COP_data_shanks, 'g');
hold on;
plot(AP_COP_complete_ts.time(peaks_APA_AP_COP), value_APA_AP_COP , 'r.');

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

% -------------Correlation between peak AP_COP and peak Acc trunk----------
figure ()

plot(value_APA_trunk_X, value_APA_AP_COP, '.r');
title('Linear correlation between peak AP_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('AP_COP (mmm)');

% -------------------------------------------------------------------------
% 2.3) Determine when th second step (when patient goes down from
% plateform) starts with left or right foot.
% -------------------------------------------------------------------------
% When the patient starts to step with left foot, last value of the cycle
% of ML COP is positive because the step finishes with the right foot, so
% the ML COP is shifted toward right (positive value). 
% When the patient starts with right foot, the last value of the ML COP
% cycle will be negative.

ML_COP_data_shanks =  reshape(ML_COP_data(3, 1, :), ...
                [1, max(size(ML_COP_data))]);

cycle_start_right = find(ML_COP_data_shanks(final_second_step_edge) < 0);
%cycle_start_left = find(ML_COP_data_shanks(final_second_step_edge) > 0);


% -------------------------------------------------------------------------
% 2.4) Find the APA peaks of ML COP and trunk acceleration of Y-axis.
% -------------------------------------------------------------------------

% If the patient starts to walk with the left foot, the APA peak is negative
% (in ML COP and Acc trunk of the Y-axis). 
% If the patient starts to walk with the right foot, the APA peak is
% positive (in ML COP and Acc trunk of the Y-axis).

% Obtain the trunk signal of Y-axis (Medio-Lateral movement). 
a_trunk_data_Y = reshape(a_trunk_data(2, 1, :), ...
                [1, max(size(a_trunk_data))]);

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point roundoff error.
    initcross_trunk = find(abs(a_trunk_complete_ts.time - initcross(k)) < 0.001);
    initcross_COP = find(abs(ML_COP_complete_ts.time - initcross(k)) < 0.001);

    finalcross_trunk = find(abs(a_trunk_complete_ts.time - finalcross(k)) < 0.001);
    finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k)) < 0.001);
    
    % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a positive peak.
        
        % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                a_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));

        % Store the index of the last negative peak.                                      
        peaks_APA_trunk_Y(k) = find(a_trunk_data_Y(...
                        initcross_trunk:finalcross_trunk)== max(...
                        neg_peak_values), 1, 'last') + initcross_trunk - 1;

         % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                ML_COP_data_shanks(...
                                initcross_COP:finalcross_COP));

        % Store the index of the last negative peak.                                      
        peaks_APA_ML_COP(k) = find(ML_COP_data_shanks(...
                        initcross_COP:finalcross_COP)== max(...
                        neg_peak_values), 1, 'last') + initcross_COP - 1;
        
    else % Look for a negative peak.
        % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -a_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));

        % Store the index of the last negative peak.                                      
        peaks_APA_trunk_Y(k) = find(a_trunk_data_Y(...
                        initcross_trunk:finalcross_trunk)== -max(...
                        neg_peak_values), 1, 'last') + initcross_trunk - 1;

         % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -ML_COP_data_shanks(...
                                initcross_COP:finalcross_COP));

        % Store the index of the last negative peak.                                      
        peaks_APA_ML_COP(k) = find(ML_COP_data_shanks(...
                        initcross_COP:finalcross_COP)== -max(...
                        neg_peak_values), 1, 'last') + initcross_COP - 1;
                    
    end               
end

% Calculate the value of the APA peaks in each cycle in the trunk 
% acceleration and the ML COP.
value_APA_trunk_Y = a_trunk_data_Y(peaks_APA_trunk_Y);
value_APA_ML_COP = ML_COP_data_shanks(peaks_APA_ML_COP);

% Calculate the linear correlation between the peaks detected with the 
% acceleration signal of the trunk and ML COP.
value_APA_trunk_Y_c = abs(value_APA_trunk_Y);
value_APA_ML_COP_c = abs( value_APA_ML_COP);
[corr_trunk_ML, prob_trunk_ML] = corr(value_APA_trunk_Y_c',...
                                    value_APA_ML_COP_c');


%------------------------------- Plots-------------------------------------

%-------------Peaks of APA in trunk acceleration and ML COP----------------
figure ();
subplot(2,1,1)
plot(a_trunk_complete_ts.time, a_trunk_data_Y, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_trunk_Y), value_APA_trunk_Y , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');
title(['Acceleration of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Acceleration (g)');


subplot(2,1,2)
plot(ML_COP_complete_ts.time, ML_COP_data_shanks, 'g');
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

% -------------Correlation between peak ML_COP and peak Acc trunk----------

figure ()

plot(value_APA_trunk_Y_c, value_APA_ML_COP_c, '.r');
title('Linear correlation between peak ML_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('ML_COP (mmm)');


% -------------------------------------------------------------------------
% 3) Detect APA in trunk signal of the gyroscope.
% -------------------------------------------------------------------------

