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

init_second_step = edges_right(3:4:length(edges_right));

force_sum_complete_ts = append(force_sum_ts{1, :});
force_sum_data = force_sum_complete_ts.data(3,1,:);
force_sum_data = reshape(force_sum_data, 1, length(force_sum_data));

final_second_step = find(diff(force_sum_data > 100)~=0);
final_second_step = final_second_step(2:2:length(final_second_step));
final_second_step = force_sum_complete_ts.time(final_second_step);

% Extract data from timeseries for plot.
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

% -------------------------------------------------------------------------
% Plots
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
% 1.1) Differences when the patient starts with left or right foot.
% -------------------------------------------------------------------------

% Determine the cycles when patient start with left or right foot.

% Beginning with the left foot.
cycle_start_left = find(sync_peaks_l < sync_peaks_r);
a_trunk_left = a_trunk(cycle_start_left);

% Inicialization of the variables
a_trunk_left_X = zeros(length(cycle_start_left),5000); 
a_trunk_left_Y = zeros(length(cycle_start_left),5000);

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
