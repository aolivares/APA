clear all; close all; clc;

% Load data.
load('data.mat');

% Compute pitch angles with EKF.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, p, X] = ...
            fusion_EKF(g_Y_right_thigh_1_C', ...
                       g_Y_right_shank_1_C', ...
                       a_X_right_thigh_1_C', ...
                       a_Z_right_thigh_1_C', ...
                       a_X_right_shank_1_C', ...
                       a_Z_right_shank_1_C', ...
                       f, 0.4, 0.3);

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
   
% matlab2tikz('../tikz/experiment_3.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');

% Plot: State vector.
n1 = 1;
n2 = length(a_X_left_shank_1_C);
figure4=figure(4);
plot(time(n1:n2), X);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('States');
l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
       '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
       '$\alpha_2$', '$\beta_1$', '$\beta_2$');
set(l,'Interpreter','Latex');
   
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');




% Plot: Acceleration based angle estimate and marker
%       signal.
n1 = 60 * f;
n2 = 80 * f;
figure5 = figure(5);
hold on;
%plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90);
%plot(time(n1:n2), p(2,n1:n2));
plot(time(n1:n2), pitch_GKF_right_thigh(n1:n2)-90, 'linewidth', 1.5);
plot(time(n1:n2), pitch_EKF_right_thigh(n1:n2), 'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Gated Kalman filter', ...
       'Extended Kalman filter');
   
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');

% Plot: Acceleration based angle estimate and marker
%       signal.
n1 = 60 * f;
n2 = 80 * f;
figure6 = figure(6);
hold on;
%plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90);
%plot(time(n1:n2), p(2,n1:n2));
plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90, 'linewidth', 1.5);
plot(time(n1:n2), pitch_EKF_right_thigh(n1:n2)+pitch_EKF_right_shank(n1:n2), 'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Gated Kalman filter', ...
       'Extended Kalman filter');
   
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');