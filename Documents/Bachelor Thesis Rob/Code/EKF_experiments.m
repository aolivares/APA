clear all; close all; clc;

% Load existing angle estimates based on the existing
% algorithms and the Qualisys motion capture system.
load('GaitWatch_data.mat');
load('Qualisys_data.mat');

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
          
% Plot: Thigh angle estimate - acceleration-based,
%       GKF, and EKF.
n1 = 1;
n2 = 80 * f;
figure4 = figure(4);
hold on;
plot(time(n1:n2), pitch_QS_right_thigh(n1:n2) ...
                  * 180 / pi - 90);
plot(time(n1:n2), pitch_acc_right_thigh(n1:n2)-90);
plot(time(n1:n2), pitch_KF_right_thigh(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2), ...
     'linewidth', 1.5);
 
xlabel('Time $t$ in s', 'interpreter', 'latex');
ylabel(['Pitch angle $\theta_1$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Reference', 'Accelerometer-based', ...
       'Kalman filter', 'Extended Kalman filter');
   
%cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz('../tikz/experiment_4.tikz', 'height', ...
           '\figureheight', 'width', '\figurewidth');

% % Plot: Thigh angle estimate - acceleration-based,
% %       GKF, and EKF.
% n1 = 62 * f + 1;
% n2 = 76 * f;
% figure5 = figure(5);
% hold on;
% plot(time(n1:n2), Q_leg_pitch(n1:n2,2)*180/pi+mean(pitch_acc_right_thigh(300))-90);
% plot(time(n1:n2), pitch_acc_right_thigh(n1:n2)-90);
% plot(time(n1:n2), pitch_GKF_right_thigh(n1:n2)-90, ...
%      time(n1:n2), pitch_EKF_right_thigh(n1:n2), ...
%      'linewidth', 1.5);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', 'Gated Kalman filter', ...
%        'Extended Kalman filter');
%    
% %cleanfigure('minimumPointsDistance', 1);
%     
% %matlab2tikz('../tikz/experiment_5.tikz', 'height', ...
% %            '\figureheight', 'width', '\figurewidth');
% 
% % Plot: Shank angle estimate - acceleration-based,
% %       GKF, and EKF.
% n1 = 1;
% n2 = 80 * f;
% figure6 = figure(6);
% hold on;
% plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
% plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90, ...
%      time(n1:n2), pitch_EKF_right_thigh(n1:n2) + ...
%      pitch_EKF_right_shank(n1:n2), 'linewidth', 1.5);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', 'Gated Kalman filter', ...
%        'Extended Kalman filter');
%    
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_6.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
% 
% % Plot: Shank angle estimate - acceleration-based,
% %       GKF, and EKF.
% n1 = 62 * f + 1;
% n2 = 76 * f;
% figure7 = figure(7);
% hold on;
% plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
% plot(time(n1:n2), pitch_GKF_right_shank(n1:n2)-90, ...
%      time(n1:n2), pitch_EKF_right_thigh(n1:n2) + ...
%      pitch_EKF_right_shank(n1:n2), 'linewidth', 1.5);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', 'Gated Kalman filter', ...
%        'Extended Kalman filter');
%    
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_7.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');       
%         
% % Plot: State vector.
% n1 = 1;
% n2 = 80 * f;
% figure8 = figure(8);
% hold on;
% plot(time(n1:n2), x(:, n1:n2));
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('States');
% l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
%          '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
%          '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
% set(l, 'Interpreter', 'Latex');
%    
% cleanfigure('minimumPointsDistance', 10);
%     
% matlab2tikz('../tikz/experiment_8.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth'); 
%         
% % Plot: State vector.
% n1 = 1;
% n2 = len;
% figure9 = figure(9);
% hold on;
% plot(time(n1:n2), x(:, n1:n2));
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('States');
% l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
%          '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
%          '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
% set(l, 'Interpreter', 'Latex');
%    
% cleanfigure('minimumPointsDistance', 10);
%     
% matlab2tikz('../tikz/experiment_9.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
%         
% % Plot: Acceleration-based pitch angle shank - corrected.
% n1 = 62 * f + 1;
% n2 = 76 * f;
% figure10 = figure(10);
% hold on;
% plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
% plot(time(n1:n2), theta12_c(n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', ...
%        'Accelerometer based - corrected');
%    
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_10.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
% 
% % Plot: Acceleration-based pitch angle shank - corrected.
% n1 = 67 * f + 1;
% n2 = 70 * f;
% figure11 = figure(11);
% hold on;
% plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
% plot(time(n1:n2), theta12_c(n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
%         '$^{\circ}$'], 'interpreter','latex');
% legend('Accelerometer-based', ...
%        'Accelerometer based - corrected');
%    
% cleanfigure('minimumPointsDistance', 1);
%     
% matlab2tikz('../tikz/experiment_11.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
%         
% % Plot: Acceleration in x-direction that sensor 2 will 
% %       see due to motion.
% n1 = 62 * f + 1;
% n2 = 76 * f;
% figure12 = figure(12);
% hold on;
% plot(time(n1:n2), a_X_right_shank_1_C(n1:n2), ...
%      time(n1:n2), a_m(1, n1:n2));
% plot(time(n1:n2), a_X_right_shank_1_C(n1:n2)' - ...
%                   a_m(1, n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('Acceleration $a_x$ in g', ...
%        'interpreter', 'latex');
% legend('Acceleration in x-direction', ...
%        'Acceleration in x-direction due to motion', ...
%        'Acceleration in x-direction - corrected');
%     
% matlab2tikz('../tikz/experiment_12.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
% 
% % Plot: Acceleration in x-direction that sensor 2 will 
% %       see due to motion.
% n1 = 67 * f + 1;
% n2 = 70 * f;
% figure13 = figure(13);
% hold on;
% plot(time(n1:n2), a_X_right_shank_1_C(n1:n2), ...
%      time(n1:n2), a_m(1, n1:n2));
% plot(time(n1:n2), a_X_right_shank_1_C(n1:n2)' - ...
%                   a_m(1, n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('Acceleration $a_x$ in g', ...
%        'interpreter', 'latex');
% legend('Acceleration in x-direction', ...
%        'Acceleration in x-direction due to motion', ...
%        'Acceleration in x-direction - corrected');
%     
% matlab2tikz('../tikz/experiment_13.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
%         
% 
% % Plot: Acceleration in z-direction that sensor 2 will 
% %       see due to motion.
% n1 = 62 * f + 1;
% n2 = 76 * f;
% figure14 = figure(14);
% hold on;
% plot(time(n1:n2), a_Z_right_shank_1_C(n1:n2), ...
%      time(n1:n2), a_m(3, n1:n2));
% plot(time(n1:n2), a_Z_right_shank_1_C(n1:n2)' - ...
%                   a_m(3, n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('Acceleration $a_z$ in g', ...
%        'interpreter', 'latex');
% legend('Acceleration in z-direction', ...
%        'Acceleration in z-direction due to motion', ...
%        'Acceleration in z-direction - corrected');
%     
% matlab2tikz('../tikz/experiment_14.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
% 
% % Plot: Acceleration in z-direction that sensor 2 will 
% %       see due to motion.
% n1 = 67 * f + 1;
% n2 = 70 * f;
% figure15 = figure(15);
% hold on;
% plot(time(n1:n2), a_Z_right_shank_1_C(n1:n2), ...
%      time(n1:n2), a_m(3, n1:n2));
% plot(time(n1:n2), a_Z_right_shank_1_C(n1:n2)' - ...
%                   a_m(3, n1:n2), 'linewidth', 1);
%  
% xlabel('Time $t$ in s', 'interpreter','latex');
% ylabel('Acceleration $a_z$ in g', ...
%        'interpreter', 'latex');
% legend('Acceleration in z-direction', ...
%        'Acceleration in z-direction due to motion', ...
%        'Acceleration in z-direction - corrected');
%     
% matlab2tikz('../tikz/experiment_15.tikz', 'height', ...
%             '\figureheight', 'width', '\figurewidth');
%         
% 
%                 