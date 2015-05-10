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
showPlotsCheck = 'no';
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


% -------------------------------------------------------------------------
% 2.2) Plots
% -------------------------------------------------------------------------

% --------------- Acceleration in shanks and force-------------------------
if strcmpi(showPlotsCheck,'yes')
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
end

% --------------------Interesting cycles (Zoom)----------------------------
% ------------------------> 6º cycle
if strcmpi(showPlotsCheck,'yes')
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
end

% -------------------------------------------------------------------------
% 3) Differences when the patient starts with left or right foot in the 
% first step to see the differences.
% -------------------------------------------------------------------------

% % Determine the cycles when patient start with left or right foot.
% 
% % Beginning the step with the left foot.
% cycle_start_left = find(sync_peaks_l < sync_peaks_r);
% a_trunk_left = a_trunk(cycle_start_left);
% 
% % Inicialization of the variables
% a_trunk_left_X = zeros(length(cycle_start_left),5000); 
% a_trunk_left_Y = zeros(length(cycle_start_left),5000);
% 
% % We group the cycles start with the left foot.
% for i = 1:length(cycle_start_left)
%     
%    a_trunk_left_data = a_trunk_left{i}.data;
%    
%    a_trunk_left_data_X = reshape(a_trunk_left_data(1, 1, :), ...
%        [1, max(size(a_trunk_left_data))]);
%    a_trunk_left_X (i,1:length(a_trunk_left_data_X)) = a_trunk_left_data_X;
%    
%    a_trunk_left_data_Y = reshape(a_trunk_left_data(2, 1, :), ...
%        [1, max(size(a_trunk_left_data))]);
%    a_trunk_left_Y (i,1:length(a_trunk_left_data_Y)) = a_trunk_left_data_Y;  
%     
% end
% 
% figure()
% plot(a_trunk_left_Y(2:3,:)');
% title('Pattern of the Acceleration signal of the trunk in the Y-axe when step happend with the left foot');
% axis([0, 2000, -0.5, 0.5]);
% 
% 
% % Beginning with the right foot.
% cycle_start_right = find(sync_peaks_r < sync_peaks_l);
% a_trunk_right = a_trunk(cycle_start_right);
% 
% 
% % Inicialization of the variables.
% a_trunk_right_X = zeros(length(cycle_start_left),5000); 
% a_trunk_right_Y = zeros(length(cycle_start_left),5000);
% 
% % We group the cycles start with the left foot.
% for i = 1:length(cycle_start_right)
%     
%    a_trunk_right_data = a_trunk_right{i}.data;
%    
%    a_trunk_right_data_X = reshape(a_trunk_right_data(1, 1, :), ...
%        [1, max(size(a_trunk_right_data))]);
%    a_trunk_right_X (i,1:length(a_trunk_right_data_X)) = a_trunk_right_data_X;
%    
%    a_trunk_right_data_Y = reshape(a_trunk_right_data(2, 1, :), ...
%        [1, max(size(a_trunk_right_data))]);
%    a_trunk_right_Y (i,1:length(a_trunk_right_data_Y)) = a_trunk_right_data_Y;  
%     
% end
% 
% figure()
% plot(a_trunk_right_Y(1:2,:)');
% title('Pattern of the Acceleration signal of the trunk in the Y-axe when step happend with the right foot');
% axis([0, 2000, -0.5, 0.5]);
% 
% figure ()
% plot(a_trunk_right_X');
% title('Pattern of the Acceleration signal of the trunk in the X-axe');
% axis([0, 2000, 0, 1]);

% -------------------------------------------------------------------------
% 4) Detect APA in trunk signal.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 4.1) Application of Activity detectors in trunk signal
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
% 4.2) Find the first 'negative' peak of each cycle in the signal 
%      of the x-axis of the acceleration of the trunk, that is, the point 
%      in time when the patient shifted backward. 
%      It is done the same with the syroscope signal.
% -------------------------------------------------------------------------

% Determine the initial and end point (time) of each interval where 
% we need to find the peaks, i.e, the second activity period of each cycle.
initcross = time_GW(init_second_step);
finalcross = final_second_step';

% Obtain the trunk signal of X-axis (Antereo-Posterior movement) and AP COP
% of both feet.
a_trunk_data_X = reshape(a_trunk_data(1, 1, :), ...
                [1, max(size(a_trunk_data))]);
g_trunk_data_X = reshape(g_trunk_data(1, 1, :), ...
                [1, max(size(g_trunk_data))]);
AP_COP_data_shanks = reshape(AP_COP_data(3, 1, :), ...
                [1, max(size(AP_COP_data))]);
            
% Obtain the samples when the patient goes up toward the platform.
init_FP_edge = edges_FP_signal(1:2:...
                        length(edges_FP_signal));
init_FP_edge = AP_COP_complete_ts.time (init_FP_edge);

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point roundoff error.
    initcross_trunk = find(abs(a_trunk_complete_ts.time - initcross(k)) ...
                        < 0.001);
    initcross_COP = find(abs(AP_COP_complete_ts.time - initcross(k))...
                        < 0.001);
    
    
    % Calculate the value of the first point of the second step to obtain 
    % afterwards the height of the COP peak.
    % The same considerating the point when the patient goes up the
    % plateform.
    initcross_trunk_complete(k) = initcross_trunk;
    init_FP = find(abs(AP_COP_complete_ts.time - init_FP_edge(k)) < 0.001);
    value_initcross_AP(k) = AP_COP_data_shanks(initcross_COP);
    value_init_FP_AP(k) = AP_COP_data_shanks(init_FP + 1);
    
    finalcross_trunk = find(abs(a_trunk_complete_ts.time - finalcross(k))...
                        < 0.001);
    finalcross_COP = find(abs(AP_COP_complete_ts.time - finalcross(k))...
                    < 0.001);
    
    % Find all peaks in each interval.
    [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -a_trunk_data_X(...
                            initcross_trunk:finalcross_trunk));

    % Store the index of the last negative peak.                                      
    peaks_APA_trunk_X(k) = find(a_trunk_data_X(...
                    initcross_trunk:finalcross_trunk)== -max(...
                    neg_peak_values), 1) + initcross_trunk - 1;
                
     % Find all peaks in each interval of the gyroscope signal.                     
     [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -g_trunk_data_X(...
                            initcross_trunk:finalcross_trunk));
                                       
    % Store the index of the last positive peak in the gyroscope signal.                                    
    peaks_APA_trunk_Gyro_X(k) = find(g_trunk_data_X(...
                    initcross_trunk:finalcross_trunk)== -max(...
                    neg_peak_values), 1) + initcross_trunk - 1;
                
     % Find all peaks in each interval of AP COP signal.
    [neg_peak_values, neg_peak_locations] = findpeaks(...
                            -AP_COP_data_shanks(...
                            initcross_COP:finalcross_COP));

    % Store the index of the first negative peak.                                      
    peaks_APA_AP_COP(k) = find(AP_COP_data_shanks(...
                    initcross_COP:finalcross_COP)== -max(...
                    neg_peak_values), 1) + initcross_COP - 1;
                
    % Find the peaks before and after the APA peak to calculate the height
    % peak.
    maximuns_peaks_first (k)= max(AP_COP_data_shanks(initcross_COP:peaks_APA_AP_COP (k) ));
    maximuns_peaks_last (k) = max(AP_COP_data_shanks(peaks_APA_AP_COP(k) :finalcross_COP));
%   
    % Segmentation of the signals with a Win= 0.25s
    fs = 200;
    Win_trunk = 0.25 * fs;
    num_segment = round((finalcross_trunk - initcross_trunk)/Win_trunk);
    Win_COP = (finalcross - initcross)/num_segment;
    
     for l=1:num_segment
         
        % Obtain every segment of the interval of trunk acc signal.
        segment_trunk  = a_trunk_data_X(initcross_trunk + (l-1)*Win_trunk:...
            initcross_trunk + l*Win_trunk);
        
        % Mean and variance of every segment
        mean_trunk (k,l) = mean(segment_trunk);
        var_trunk (k,l) = var (segment_trunk);
        
        % Obtain every segment of the interval of COP AP signal.
        segment_COP_AP = AP_COP_data_shanks(round(initcross_COP + (l-1)*Win_COP):...
                    round(initcross_COP+ l*Win_COP));
                
        % Mean and variance of every segment.
        mean_COP_AP(k,l) = mean(segment_COP_AP);
        var_COP_AP (k,l) = var (segment_COP_AP); 
   
     end
    
end

% Calculate the value of the APA peaks in each cycle in the trunk 
% acceleration and the AP COP. 
value_APA_trunk_X = a_trunk_data_X(peaks_APA_trunk_X);
value_APA_trunk_Gyro_X = g_trunk_data_X(peaks_APA_trunk_Gyro_X);
value_APA_AP_COP = AP_COP_data_shanks(peaks_APA_AP_COP);
% 
% % Calculate the linear correlation between the peaks detected with the 
% % acceleration signal of the trunk and AP COP.
% [corr_trunk_AP, prob_trunk_AP] = corr(value_APA_trunk_X', value_APA_AP_COP');
% 
% % Calculate the correlation between the height of the AP COP peaks.
value_APA_AP_COP_1 = abs(value_APA_AP_COP - value_initcross_AP);
% [corr_trunk_AP_2, prob_trunk_AP_2] = corr(value_APA_trunk_X',...
%                                     value_APA_AP_COP_2');
%                                 
% value_APA_AP_COP_3 = value_APA_AP_COP - value_init_FP_AP;
% [corr_trunk_AP_3, prob_trunk_AP_3] = corr(value_APA_trunk_X',...
%                                     value_APA_AP_COP_3');

% This is the most iteresting correlation.                               
value_APA_trunk_X_1 = abs(value_APA_trunk_X - a_trunk_data_X(initcross_trunk_complete));
[corr_trunk_AP_1, prob_trunk_AP_1] = corr(value_APA_trunk_X_1',...
                                    value_APA_AP_COP_1');

value_APA_trunk_X_2 = abs(value_APA_trunk_X - mode(a_trunk_data_X));
value_APA_AP_COP_2 = abs(value_APA_AP_COP - maximuns_peaks_first);
[corr_trunk_AP_2, prob_trunk_AP_2] = corr(value_APA_trunk_X_2',...
                                    value_APA_AP_COP_2');
                                
% Calculate the correlation between the Gyro trunk signal peaks and AP COP.
% peaks.
value_APA_trunk_Gyro_X_1 = abs(value_APA_trunk_Gyro_X - g_trunk_data_X(initcross_trunk_complete));
[corr_trunk_AP_Gyro_1, prob_trunk_AP_Gyro_1] = corr(value_APA_trunk_Gyro_X_1',...
                                            value_APA_AP_COP_1');
                                        
value_APA_trunk_Gyro_X_2 = abs(value_APA_trunk_Gyro_X - mode(g_trunk_data_X));                                        
[corr_trunk_AP_Gyro_2, prob_trunk_AP_Gyro_2] = corr(value_APA_trunk_Gyro_X_2',...
                                            value_APA_AP_COP_2');
                                        
% Calculate the linear correlation between the patter that characterise the
% patient. This correlation shows the similitudes between the trunk and
% AP_COP signals.

mean_mean_trunk = mean(mean_trunk);
mean_mean_COP_AP = mean(mean_COP_AP);
[corr_mean, prob_mean] = corr(mean_trunk', mean_COP_AP');
corr_mean_mean= mean(corr_mean(:,1));
prob_mean_mean = mean(prob_mean(:,1));
                                
%--------------------------------------------------------------------------
% Plots
%--------------------------------------------------------------------------

%-------------Peaks of APA in trunk acceleration and AP COP----------------
if strcmpi(showPlotsAPA,'yes')
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
end

% -------------Correlation between peak AP_COP and peak Acc trunk----------
if strcmpi(showPlotsCorr,'yes')
figure ()

plot(value_APA_trunk_X_1, value_APA_AP_COP_2, '.r');
title('Linear correlation between peak AP_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('AP_COP (mmm)');
end


% -------------------------------------------------------------------------
% 4.3) Determine when th second step (when patient goes down from
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

% -------------------------------------------------------------------------
% 4.4) Find the APA peaks of ML COP and trunk acceleration of Y-axis.
% -------------------------------------------------------------------------

% If the patient starts to walk with the left foot, the APA peak is negative
% (in ML COP and Acc trunk of the Y-axis). 
% If the patient starts to walk with the right foot, the APA peak is
% positive (in ML COP and Acc trunk of the Y-axis).

% Obtain the trunk signal of Y-axis (Medio-Lateral movement). 
a_trunk_data_Y = reshape(a_trunk_data(2, 1, :), ...
                [1, max(size(a_trunk_data))]);
g_trunk_data_Y = reshape(g_trunk_data(2, 1, :), ...
                [1, max(size(g_trunk_data))]);
            
% Calculate the peaks in each interval.
for k = 1:length(initcross)
    
    % To find a noninteger value, we use a tolerance value based on our data.
    % Otherwise, the result is sometimes an empty matrix due to 
    % floating-point roundoff error.
    initcross_trunk = find(abs(a_trunk_complete_ts.time - initcross(k))...
                        < 0.001);
    initcross_COP = find(abs(ML_COP_complete_ts.time - initcross(k))...
                        < 0.001);
    
    % Calculate the value of the first point to obtain afterwards the
    % height of the COP peak.
    value_initcross_ML(k) = ML_COP_data_shanks(initcross_COP);

    finalcross_trunk = find(abs(a_trunk_complete_ts.time - finalcross(k))...
                        < 0.001);
    finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k))...
                        < 0.001);
    
    % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a positive peak.
        
        % Find all peaks in each interval before a negative peak.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -a_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));
                                   
        peaks_APA_trunk_Y_neg = find(a_trunk_data_Y(...
                        initcross_trunk:finalcross_trunk)== -max(...
                        neg_peak_values), 1) + initcross_trunk - 1;
                    
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                a_trunk_data_Y(...
                                initcross_trunk:peaks_APA_trunk_Y_neg));

        % Store the index of the last negative peak.                                      
        peaks_APA_trunk_Y(k) = find(a_trunk_data_Y(...
                        initcross_trunk:finalcross_trunk)== max(...
                        pos_peak_values), 1) + initcross_trunk - 1;

         % Find all peaks in each interval.
         
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                ML_COP_data_shanks(...
                                initcross_COP:finalcross_COP));

        % Store the index of the last positive peak.                                      
        peaks_APA_ML_COP(k) = find(ML_COP_data_shanks(...
                        initcross_COP:finalcross_COP)== max(...
                        pos_peak_values), 1) + initcross_COP - 1;
                    
        % Find positives peaks in gyroscope trunk signal before the largest 
        % negative peak.
         [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -g_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));
                            
         peaks_APA_trunk_Gyro_Y_neg = find(g_trunk_data_Y(...
                    initcross_trunk:finalcross_trunk)== -max(...
                    neg_peak_values), 1) + initcross_trunk - 1;     
    
         [pos_peak_values, pos_peak_locations] = findpeaks(...
                                g_trunk_data_Y(...
                                initcross_trunk:peaks_APA_trunk_Gyro_Y_neg));
                            
         % Check if there are positive peaks detected.
         if isnan(pos_peak_values)
            [pos_peak_values, pos_peak_locations] = findpeaks(...
                                g_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));
         end
                            
       % Store the index of the last positive peak.                                      
        peaks_APA_trunk_Gyro_Y(k) = find(g_trunk_data_Y(...
                    initcross_trunk:finalcross_trunk)== max(...
                    pos_peak_values), 1) + initcross_trunk - 1;
        
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
                        neg_peak_values), 1) + initcross_COP - 1;
                    
        % Find negatives peaks in gyroscope trunk signal before the largest
        % positives peaks.
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                g_trunk_data_Y(...
                                initcross_trunk:finalcross_trunk));
                                           
        peaks_APA_trunk_Gyro_Y_pos = find(g_trunk_data_Y(...
                    initcross_trunk:finalcross_trunk)== max(...
                    pos_peak_values), 1) + initcross_trunk - 1;
                
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                -g_trunk_data_Y(...
                                initcross_trunk:peaks_APA_trunk_Gyro_Y_pos));
                            
         % Check if there are negative peaks detected.
         if isnan(neg_peak_values)else
         [neg_peak_values, neg_peak_locations] = findpeaks(...
                    -g_trunk_data_Y(...
                    initcross_trunk:finalcross_trunk)); 
         end

        % Store the index of the last negative peak.                                      
         peaks_APA_trunk_Gyro_Y(k) = find(g_trunk_data_Y(...
                    initcross_trunk:finalcross_trunk)== -max(...
                    neg_peak_values), 1) + initcross_trunk - 1;
                    
    end
                    
end

% Calculate the value of the APA peaks in each cycle in the trunk 
% acceleration and the ML COP.
value_APA_trunk_Y = a_trunk_data_Y(peaks_APA_trunk_Y);
value_APA_trunk_Gyro_Y = g_trunk_data_Y(peaks_APA_trunk_Gyro_Y);
value_APA_ML_COP = ML_COP_data_shanks(peaks_APA_ML_COP);

% Calculate the linear correlation between the peaks detected with the 
% acceleration signal of the trunk and ML COP.
value_APA_trunk_Y_c = abs(value_APA_trunk_Y);
value_APA_ML_COP_c = abs( value_APA_ML_COP);
[corr_trunk_ML, prob_trunk_ML] = corr(value_APA_trunk_Y_c',...
                                    value_APA_ML_COP_c');

                                
% Calculate the correlation between the height of the ML COP peaks.
value_APA_ML_COP_c_2 = value_APA_ML_COP_c - abs(value_initcross_ML);
[corr_trunk_ML_2, prob_trunk_ML_2] = corr(value_APA_trunk_Y_c',...
                                    value_APA_ML_COP_c_2');
                                
value_APA_trunk_Y_1 = abs(value_APA_trunk_Y_c - a_trunk_data_Y(...
                        initcross_trunk_complete));
[corr_trunk_ML_3, prob_trunk_ML_3] = corr(value_APA_trunk_Y_1',...
                                    value_APA_ML_COP_c_2');
                                
value_APA_trunk_Y_2 = abs(value_APA_trunk_Y_c - mode(a_trunk_data_Y));
[corr_trunk_ML_4, prob_trunk_ML_4] = corr(value_APA_trunk_Y_2',...
                                    value_APA_ML_COP_c');
% Correlation gyroscope signal.
value_APA_trunk_Gyro_Y_1 = abs(value_APA_trunk_Gyro_Y - g_trunk_data_Y(...
                            initcross_trunk_complete));
[corr_trunk_ML_Gyro, prob_trunk_ML_Gyro] = corr(abs(value_APA_trunk_Gyro_Y_1)',...
                                    value_APA_ML_COP_c_2');
                                

%------------------------------- Plots-------------------------------------

%-------------Peaks of APA in trunk acceleration and ML COP----------------
if strcmpi(showPlotsAPA,'yes')
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
end

% -------------Correlation between peak ML_COP and peak Acc trunk----------
if strcmpi(showPlotsCorr,'yes')
figure ()

plot(value_APA_trunk_Y_1,value_APA_ML_COP_c_2, '.r');
title('Linear correlation between peak ML_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('ML_COP (mmm)');

end

% --------------------Angular Velocity in trunk and COP--------------------
if strcmpi(showPlotsCheck,'yes')
figure()
subplot(4, 1, 1);
plot(g_trunk_complete_ts.time, g_trunk_data_X )

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Angular Velocity of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');


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
plot(g_trunk_complete_ts.time, g_trunk_data_Y);

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Angular velocity of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time']);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');

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
end
% peaks--------------------------------------------------------------------
if strcmpi(showPlotsAPA,'yes')
figure ();
subplot(2,1,1)
plot(g_trunk_complete_ts.time, g_trunk_data_X, 'g');
hold on;
plot(g_trunk_complete_ts.time(peaks_APA_trunk_Gyro_X),value_APA_trunk_Gyro_X, 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');
title(['Angular Velocity of the X-axis of the trunk with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');


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

figure ();
subplot(2,1,1)
plot(g_trunk_complete_ts.time, g_trunk_data_Y, 'g');
hold on;
plot(g_trunk_complete_ts.time(peaks_APA_trunk_Gyro_Y),value_APA_trunk_Gyro_Y, 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');
title(['Angular Velocity of the Y-axis of the trunk with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');


subplot(2,1,2)
plot(ML_COP_complete_ts.time, ML_COP_data_shanks, 'g');
hold on;
plot(AP_COP_complete_ts.time(peaks_APA_ML_COP), value_APA_ML_COP , 'r.');

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


