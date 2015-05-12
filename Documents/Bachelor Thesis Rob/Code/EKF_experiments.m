clear all; close all; clc;

% Load existing angle estimates based on the existing
% algorithms and the Qualisys motion capture system.
load('GaitWatch_data_2.mat');
load('Qualisys_data_2.mat');

% Initialise number of figure.
n = 4;

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
%       KF, and EKF.
n1 = 1;
n2 = 24 * f;
figure(n);
hold on;
plot(time(n1:n2), pitch_QS_right_thigh(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2), pitch_acc_right_thigh(n1:n2) - 90);
plot(time(n1:n2), pitch_KF_right_thigh(n1:n2) - 90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2), ...
     'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter', 'latex');
ylabel(['Pitch angle $\theta_1$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Reference', 'Accelerometer-based', ...
       'Kalman filter', 'Extended Kalman filter');
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
       
n = n + 1;

% Plot: Shank angle estimate - acceleration-based,
%       KF, and EKF.
n1 = 1;
n2 = 24 * f;
figure(n);
hold on;
plot(time(n1:n2), pitch_QS_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2), pitch_acc_right_shank(n1:n2)-90);
plot(time(n1:n2), pitch_KF_right_shank(n1:n2)-90, ...
     time(n1:n2), pitch_EKF_right_thigh(n1:n2) ...
     + pitch_EKF_right_shank(n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter', 'latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter','latex');
legend('Reference', 'Accelerometer-based', ...
       'Kalman filter', 'Extended Kalman filter');
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
        
n = n + 1;        

% Compute root-mean-square error.       
RMSE_acc = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
                - pitch_acc_right_thigh(n1:n2)').^2));       
RMSE_KF = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
               - pitch_KF_right_thigh(n1:n2)').^2));
RMSE_EKF = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
             - 90 - pitch_EKF_right_thigh(n1:n2)).^2));
RMSE = [RMSE_acc, RMSE_KF, RMSE_EKF];

RMSE_acc = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
                - pitch_acc_right_shank(n1:n2)').^2));       
RMSE_KF = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
               - pitch_KF_right_shank(n1:n2)').^2));
RMSE_EKF = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
             - 90 - pitch_EKF_right_thigh(n1:n2) ...
             - pitch_EKF_right_shank(n1:n2)).^2));
RMSE = [RMSE; RMSE_acc, RMSE_KF, RMSE_EKF];

figure(n);
b = bar(RMSE, 0.3);
offset = 0.5;
yb = cat(1, b.YData); 
xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
hold on;
for i = 1:2  
   for j = 1:3
        text(xb(j, i),yb(j, i)+offset, num2str(...
        RMSE(i, j),'$%0.2f$'), 'rotation', 90, ...
        'interpreter','latex');
   end
end

b(1).FaceColor = [0.8500    0.3250    0.0980];
b(2).FaceColor = [0.9290    0.6940    0.1250];
b(3).FaceColor = [0.4940    0.1840    0.5560];

text(1,20, ['$RMSE_{thigh} + RMSE_{shank} = ', ...
            num2str(RMSE(3)+RMSE(6),'%0.2f$')], ...
            'interpreter','latex');
        
        RMSE(3) + RMSE(6)

ylim([0, max(max(RMSE)) + 6]);
ylabel('Root-mean-square error in $^{\circ}$', ...
       'interpreter','latex');

labels = {'Thigh', 'Shank'};
format_ticks(gca, labels, [], [], [], 0);

legend('Acceleration-based', 'Kalman filter', ...
          'Extended Kalman filter');

matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
       
n = n + 1;

% Plot: State vector.
n1 = 1;
n2 = 20 * f;
figure(n);
hold on;
plot(time(n1:n2), x(:, n1:n2));
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('States');
l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
         '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
         '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
set(l, 'Interpreter', 'Latex');
   
cleanfigure('minimumPointsDistance', 0.5);
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth'); 

n = n + 1;        
        
%Plot: State vector.
n1 = 1;
n2 = 40 * f;
figure(n);
hold on;
plot(time(n1:n2), x(:, n1:n2));
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('States');
l=legend('$x$', '$z$', '$\theta_1$', '$\omega_1$', ...
         '$\alpha_1$', '$\theta_2$', '$\omega_2$', ...
         '$\alpha_2$', '$\beta_1$',  '$\beta_2$');
set(l, 'Interpreter', 'Latex');
   
cleanfigure('minimumPointsDistance', 1);
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
         
n = n + 1;
         
% Plot: Acceleration-based pitch angle shank - corrected.
n1 = 4 * f + 1;
n2 = 20 * f;
figure(n);
hold on;
plot(time(n1:n2) - 4, pitch_QS_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2) - 4, pitch_acc_right_shank(n1:n2) - 90);
plot(time(n1:n2) - 4, theta12_c(n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter', 'latex');
legend('Reference', 'Accelerometer-based', ...
       'Accelerometer based - corrected');
   
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
         
n = n + 1;

% Compute root-mean-square error.       
RMSE_acc = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
                - pitch_acc_right_shank(n1:n2)').^2));       
RMSE_acc_corr = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
                - 90 - theta12_c(n1:n2)).^2));
RMSE = [RMSE_acc, RMSE_acc_corr];

figure(n);
bar(RMSE, 0.3);
ylim([0, max(RMSE)+2]);
ylabel('Root-mean-square error in $^{\circ}$', ...
       'interpreter','latex');
labels = {'Accelerometer-based', ...
          'Accelerometer based - corrected'};

format_ticks(gca, labels, [], [], [], 0);

text(1:2, RMSE' + 0.7, num2str(RMSE','$%0.2f$'),... 
'HorizontalAlignment','center', 'interpreter','latex');

matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
       
n = n + 1;

% Plot: Acceleration in x-direction that sensor 2 will 
%       see due to motion.
n1 = 4 * f + 1;
n2 = 20 * f;
figure(n);
hold on;
plot(time(n1:n2) - 4, a_X_right_shank_1_C(n1:n2));
plot(time(n1:n2) - 4, a_m(1, n1:n2), 'linewidth', 1);
plot(time(n1:n2) - 4, a_X_right_shank_1_C(n1:n2)' - ...
                      a_m(1, n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('Acceleration $a_x$ in g', ...
       'interpreter', 'latex');
legend('Acceleration in x-direction', ...
       'Acceleration in x-direction due to motion', ...
       'Acceleration in x-direction - corrected');
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');

n = n + 1;

% Plot: Acceleration in z-direction that sensor 2 will 
%       see due to motion.
n1 = 4 * f + 1;
n2 = 20 * f;
figure(n);
hold on;
plot(time(n1:n2) - 4, a_Z_right_shank_1_C(n1:n2));
plot(time(n1:n2) - 4, a_m(3, n1:n2), 'linewidth', 1);
plot(time(n1:n2) - 4, a_Z_right_shank_1_C(n1:n2)' - ...
                      a_m(3, n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel('Acceleration $a_z$ in g', ...
       'interpreter', 'latex');
legend('Acceleration in z-direction', ...
       'Acceleration in z-direction due to motion', ...
       'Acceleration in z-direction - corrected');
    
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');

n = n + 1;
        
        