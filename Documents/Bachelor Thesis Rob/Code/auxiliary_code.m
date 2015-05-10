clear all; close all; clc;

% Load data.
load('data.mat');

% Compute length of signal vectors.
len = length(a_X_right_thigh_1_C);

% Compute pitch angles with EKF.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, p, X] = ...
            fusion_EKF(g_Y_right_thigh_1_C', ...
                       g_Y_right_shank_1_C', ...
                       a_X_right_thigh_1_C', ...
                       a_Z_right_thigh_1_C', ...
                       a_X_right_shank_1_C', ...
                       a_Z_right_shank_1_C', ...
                       f, 0.3, 0.25);

% % Plot: Acceleration-based estimate vs. integration of 
% %       the angular rate.
% n = 18 * f;
% figure1 = figure(1);
% hold on
% plot(time(1:n), pitch_acc_right_shank(1:n) - 90, ...
%      time(1:n), pitch_gyro_right_shank(1:n) - 90);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', ...
%        'Integration of angular rate');
%   
% matlab2tikz('../tikz/experiment_1.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        
% % Plot: Acceleration-based estimate vs. classic Kalman
% %       filter.
% n = 18 * f;
% figure2 = figure(2);
% hold on
% plot(time(1:n), pitch_acc_right_shank(1:n) - 90);
% plot(time(1:n), pitch_KF_right_shank(1:n) - 90, ...
%      'linewidth', 1.5);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', 'Kalman filter');
%   
% matlab2tikz('../tikz/experiment_2.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        
% % Plot: Acceleration-based estimate vs. classic Kalman
% %       filter and gated Kalman filter.
% n = 18 * f;
% figure3 = figure(3);
% hold on
% plot(time(1:n), pitch_acc_right_shank(1:n) - 90);
% plot(time(1:n), pitch_KF_right_shank(1:n) - 90);
% plot(time(1:n), pitch_GKF_right_shank(1:n) - 90, ...
%      'linewidth', 1.5);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', 'Kalman filter', ...
%        'Gated Kalman filter');
%   
% matlab2tikz('../tikz/experiment_3.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');



   