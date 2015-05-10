clear all; close all; clc;

% Load data.
load('data.mat');

% Compute length of signal vectors.
len = length(a_X_right_thigh_1_C);

% Compute pitch angles with EKF and store state vector
% at each time step in x.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, x, ...
 p] = ...
   fusion_EKF(g_Y_right_thigh_1_C', ...
              g_Y_right_shank_1_C', ...
              a_X_right_thigh_1_C', ...
              a_Z_right_thigh_1_C', ...
              a_X_right_shank_1_C', ...
              a_Z_right_shank_1_C', ...
              f, 0.3, 0.25);
          
% Plot: Thigh angle estimate - acceleration-based,
%       GKF, and EKF.
n1 = 1;
n2 = 80 * f;
figure4 = figure(4);
hold on;
plot(time(n1:n2), pitch_acc_right_thigh(n1:n2)-90);
plot(time(n1:n2), pitch_GKF_right_thigh(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2), ...
     'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Gated Kalman filter', ...
       'Extended Kalman filter');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_4.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');

% Plot: Thigh angle estimate - acceleration-based,
%       GKF, and EKF.
n1 = 62 * f + 1;
n2 = 76 * f;
figure5 = figure(5);
hold on;
plot(time(n1:n2), pitch_acc_right_thigh(n1:n2)-90);
plot(time(n1:n2), pitch_GKF_right_thigh(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2), ...
     'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Gated Kalman filter', ...
       'Extended Kalman filter');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');

% Plot: Shank angle estimate - acceleration-based,
%       GKF, and EKF.
n1 = 1;
n2 = 80 * f;
figure6 = figure(6);
hold on;
plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2) + ...
     pitch_EKF_right_shank(n1:n2), 'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Gated Kalman filter', ...
       'Extended Kalman filter');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_6.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');

% Plot: Thigh angle estimate - acceleration-based,
%       GKF, and EKF.
n1 = 62 * f + 1;
n2 = 76 * f;
figure7 = figure(7);
hold on;
plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2) + ...
     pitch_EKF_right_shank(n1:n2), 'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Gated Kalman filter', ...
       'Extended Kalman filter');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_7.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');       
        
% Plot: State vector.
n1 = 1;
n2 = 80 * f;
figure8 = figure(8);
hold on;
plot(time(n1:n2), x(:, n1:n2));
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('States');
l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
         '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
         '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
set(l, 'Interpreter', 'Latex');
   
cleanfigure('minimumPointsDistance', 10);
    
matlab2tikz('../tikz/experiment_8.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth'); 
        
% Plot: State vector.
n1 = 1;
n2 = len;
figure9 = figure(9);
hold on;
plot(time(n1:n2), x(:, n1:n2));
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('States');
l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
         '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
         '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
set(l, 'Interpreter', 'Latex');
   
cleanfigure('minimumPointsDistance', 10);
    
matlab2tikz('../tikz/experiment_9.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');
        
% Plot: Acceleration correction.
n1 = 62 * f + 1;
n2 = 76 * f;
figure10 = figure(10);
hold on;
plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
plot(time(n1:n2), p(2, n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Accelerometer based - corrected');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_10.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');

% Plot: Acceleration correction.
n1 = 67 * f + 1;
n2 = 70 * f;
figure11 = figure(11);
hold on;
plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
plot(time(n1:n2), p(2, n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Accelerometer-based', 'Accelerometer based - corrected');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_11.tikz', 'height', ...
            '\figureheight', 'width', '\figurewidth');