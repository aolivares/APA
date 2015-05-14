clear all; close all; clc;

% Load existing angle estimates based on the existing
% algorithms and the Qualisys motion capture system.
load('GaitWatch_data_2.mat');
load('Qualisys_data_2.mat');

% Compute length of the signal vectors.
len = length(a_X_right_thigh_1_C);

% Compute pitch angles with extended Kalman filter.
% Additionally, store the internal state vector at each
% time step in x, the motion based acceleration in a_m,
% and the angle estimate theta_1 + theta_2, based on 
% the corrected acceleration signal in theta12_c.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, ...
 theta12_c, a_m, x] = fusion_EKF( ...
                            g_Y_right_thigh_1_C', ...
                            g_Y_right_shank_1_C', ...
                            a_X_right_thigh_1_C', ...
                            a_Z_right_thigh_1_C', ...
                            a_X_right_shank_1_C', ...
                            a_Z_right_shank_1_C', ...
                            f, 0.3, 0.25);

% Plot: Acceleration-based estimate vs. integration of 
%       the angular rate.
n1 = 1;
n2 = 18 * f;
figure1 = figure(1);
hold on
plot(time(n1:n2), pitch_QS_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2), pitch_acc_right_shank(n1:n2) - 90);
plot(time(n1:n2), pitch_gyro_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angles $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter', 'latex');
legend('Reference', 'Accelerometer-based', ...
       'Integration of angular rate');
  
matlab2tikz('../tikz/experiment_1.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');
  
% Plot: Acceleration-based estimate vs. classic Kalman
%       filter.
n1 = 1;
n2 = 18 * f;
figure2 = figure(2);
hold on
plot(time(n1:n2), pitch_QS_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2), pitch_acc_right_shank(n1:n2) - 90);
plot(time(n1:n2), pitch_KF_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angles $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter', 'latex');
legend('Reference', 'Accelerometer-based', ...
       'Kalman filter');
  
matlab2tikz('../tikz/experiment_2.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');
        
% % Plot: Acceleration-based estimate vs. classic Kalman
% %       filter.
% n1 = 1;
% n2 = 18 * f;
% figure3 = figure(3);
% hold on
% plot(time(n1:n2), pitch_QS_right_shank(n1:n2) - 90, ...
%      'linewidth', 1);
% plot(time(n1:n2), pitch_acc_right_shank(n1:n2) - 90);
% plot(time(n1:n2), pitch_GKF_right_shank(n1:n2) - 90, ...
%      'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angles $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter', 'latex');
% legend('Reference', 'Accelerometer-based', ...
%        'Gated Kalman filter');
%   
% matlab2tikz('../tikz/experiment_3.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        

        