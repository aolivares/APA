clear all; close all; clc;

% Load data.
load('data.mat');

% Compute pitch angles with EKF.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, p] = ...
            fusion_EKF(g_Y_right_thigh_1_C', ...
                       g_Y_right_shank_1_C', ...
                       a_X_right_thigh_1_C', ...
                       a_Z_right_thigh_1_C', ...
                       a_X_right_shank_1_C', ...
                       a_Z_right_shank_1_C', ...
                       f, 0.25, 0.25);

% Plot: Acceleration-based estimate vs. integration of 
%       the angular rate.
n = 18 * f;
figure1 = figure(1);
hold on
plot(time(1:n), pitch_acc_right_shank(1:n) - 90, ...
     time(1:n), pitch_gyro_right_shank(1:n) - 90);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', ...
       'Integration of angular rate');
   
% matlab2tikz('../tikz/experiment_1.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        
% Plot: Acceleration-based estimate vs. classic Kalman
%       filter.
n = 18 * f;
figure2 = figure(2);
hold on
plot(time(1:n), pitch_acc_right_shank(1:n) - 90);
plot(time(1:n), pitch_KF_right_shank(1:n) - 90, ...
     'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Kalman filter');
   
% matlab2tikz('../tikz/experiment_2.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        
% Plot: Acceleration-based estimate vs. classic Kalman
%       filter and gated Kalman filter.
n = 18 * f;
figure3 = figure(3);
hold on
plot(time(1:n), pitch_acc_right_shank(1:n) - 90);
plot(time(1:n), pitch_KF_right_shank(1:n) - 90);
plot(time(1:n), pitch_GKF_right_shank(1:n) - 90, ...
     'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Kalman filter', ...
       'Gated Kalman filter');
   
% matlab2tikz('../tikz/experiment_3.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
        
% Plot: Acceleration-based estimate vs. classic Kalman
%       filter, gated Kalman filter, and extended
%       Kalman filter.
n = 18 * f;
figure4 = figure(4);
hold on
plot(time(1:n), pitch_acc_right_shank(1:n) - 90, ...
     time(1:n), pitch_KF_right_shank(1:n) - 90, ...
     time(1:n), pitch_GKF_right_shank(1:n) - 90);
% plot(time(1:n), p(1:n));
plot(time(1:n), pitch_EKF_right_thigh(1:n) ...
                + pitch_EKF_right_shank(1:n), ...
                'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Kalman filter', ...
       'Gated Kalman filter', 'Extended Kalman filter');
    
% matlab2tikz('../tikz/experiment_4.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');

% Plot: Acceleration-based estimate vs. classic Kalman
%       filter, gated Kalman filter, and extended
%       Kalman filter.
n = length(pitch_acc_left_thigh);
figure5 = figure(5);
hold on
plot(time(1:n), p(1:n), ...
     time(1:n), pitch_GKF_right_shank(1:n) - 90);
% plot(time(1:n), p(1:n));
plot(time(1:n), pitch_EKF_right_thigh(1:n) ...
                + pitch_EKF_right_shank(1:n), ...
                'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Gated Kalman filter', ...
       'Extended Kalman filter');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');