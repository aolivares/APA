%% 0) Initialisation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -----------------------------------------------------

clear all; close all; clc;

% Generate Tikz-pictures: Yes (1), or no (0).
tikz = 1;

% Load existing signals of the GaitWatch and the 
% Qualisys motion capture system.
load('GaitWatch_data_1.mat');
load('Qualisys_data_1.mat');
load('pmin_1.mat');

% Import GaitWatch and WaGyroMag functions library.
% All existing functions have to be called using
% either 'gw.functionName' or 'wag.functionName'.
gw = gwLibrary;
wag = wagLibrary;

% Initialise number of figure.
n = 4;

% Set value of the magnitude of the gravity vector. 
g = 9.81;

% Set the first and the last sample, that is, the 
% interval of the signal that is used for the following
% computations.
n1 = 1;
n2 = 24 * f;

%% 1) Optimise filter parameters \\\\\\\\\\\\\\\\\\\\\\
% -----------------------------------------------------

% Set RMSE offset.
rmse_offset = 1;

% ---KALMAN FILTER THIGH-------------------------------

% Set initial value of parameters;
p0_KF = [1000 0.001];
    
% Call the optimization routine.
[xmin, fmin, ct] = gw.optimize_KF( ...
            g_Y_right_thigh_1_C(n1:n2)', ...
            pitch_acc_right_thigh(n1:n2)', ...
            f, var(pitch_acc_right_thigh(n1:n2)), ...
            var(pitch_acc_right_thigh(n1:n2)), ...
            var(g_Y_right_thigh_1_C(n1:n2)), ...
            pitch_acc_right_thigh(1), ...
            pitch_QS_right_thigh(n1:n2), ...
            p0_KF, rmse_offset);

fprintf('-------------KF OPTIMISATION------------\n');
fprintf(['The optimisation process finished in %d ', ...
         'iterations.\n'], ct);
fprintf('The minimum RMSE found is: %0.4f\n', fmin);
fprintf(['Optimal parameters are: \n -Alpha: %0.4f', ...
        '\n -Beta: %0.4f\n'], xmin(1), xmin(2))
fprintf('----------------------------------------\n')

% Extract optimal parameters.
opt_alpha_KF = xmin(1);
opt_beta_KF = xmin(2);

% Compute thigh angle. 
pitch_KF_right_thigh = gw.fusion_KF( ...
        g_Y_right_thigh_1_C, pitch_acc_right_thigh, ...
        f, var(pitch_acc_right_thigh), ...
        var(pitch_acc_right_thigh),...
        var(g_Y_right_thigh_1_C), opt_alpha_KF, ...
        opt_beta_KF, pitch_acc_right_thigh(1));

% ---KALMAN FILTER SHANK-------------------------------

% Set initial value of parameters;
p0_KF = [1000 0.001];

% Call the optimization routine.
[xmin, fmin, ct] = gw.optimize_KF( ...
                g_Y_right_shank_1_C(n1:n2)', ...
                pitch_acc_right_shank(n1:n2)', f, ...
                var(pitch_acc_right_shank(n1:n2)), ...
                var(pitch_acc_right_shank(n1:n2)),...
                var(g_Y_right_shank_1_C(n1:n2)), ...
                pitch_acc_right_shank(1), ...
                pitch_QS_right_shank(n1:n2), ...
                p0_KF, rmse_offset);

fprintf('-------------KF OPTIMISATION------------\n');
fprintf(['The optimisation process finished in %d ', ...
         'iterations.\n'], ct);
fprintf('The minimum RMSE found is: %0.4f\n', fmin);
fprintf(['Optimal parameters are: \n -Alpha: %0.4f', ...
        '\n -Beta: %0.4f\n'], xmin(1), xmin(2))
fprintf('----------------------------------------\n')

% Extract optimal parameters.
opt_alpha_KF = xmin(1);
opt_beta_KF = xmin(2);
    
% Compute shank angle. 
pitch_KF_right_shank = gw.fusion_KF( ...
        g_Y_right_shank_1_C, pitch_acc_right_shank, ...
        f, var(pitch_acc_right_shank), ...
        var(pitch_acc_right_shank),...
        var(g_Y_right_shank_1_C), opt_alpha_KF, ...
        opt_beta_KF, pitch_acc_right_shank(1));

% ---EXTENDED KALMAN FILTER----------------------------

% Set initial value of parameters.
p0 = [3.5, 30, 30, 300, 0.001, 0.05];

 % Call the optimization routine.
[pmin, fmin, ct] = gw.optimize_EKF( ...
                    g_Y_right_thigh_1_C(n1:n2)', ...
                    g_Y_right_shank_1_C(n1:n2)', ...
                    a_X_right_thigh_1_C(n1:n2)', ...
                    a_Z_right_thigh_1_C(n1:n2)', ...
                    a_X_right_shank_1_C(n1:n2)', ...
                    a_Z_right_shank_1_C(n1:n2)', ...
                    f, 0.35, 0.25, ...
                    [pitch_QS_right_thigh(n1:n2); ...
                     pitch_QS_right_shank(n1:n2)] - 90, ...
                    p0, rmse_offset);
                
fprintf('-------------EKF OPTIMISATION------------\n');
fprintf(['The optimisation process finished in %d ', ...
         'iterations.\n'], ct);
fprintf('The minimum RMSE found is: %0.4f\n', fmin);
fprintf(['Optimal parameters are: \n -sigma_t1: ', ...
         '%0.4f\n -sigma_t2: %0.4f\n -sigma_b: ', ...
         '%0.4f\n -sigma_s_1: %0.4f\n -sigma_s_2:', ...
         ' %0.4f'], pmin);
fprintf('----------------------------------------\n')
%%
% Compute pitch angles with extended Kalman filter.
% Additionally, store the internal state vector at each
% time step in x, the motion based acceleration in a_m,
% and the angle estimate theta_1 + theta_2, based on 
% the corrected acceleration signal in theta12_c.
[pitch_EKF_right_thigh, pitch_EKF_right_shank, ...
 theta12_c, a_m, x, par] = fusion_EKF( ...
                            g_Y_right_thigh_1_C', ...
                            g_Y_right_shank_1_C', ...
                            a_X_right_thigh_1_C', ...
                            a_Z_right_thigh_1_C', ...
                            a_X_right_shank_1_C', ...
                            a_Z_right_shank_1_C', ...
                            f, 0.35, 0.25, pmin);

save ('Data_1/pmin_1', 'pmin')  
save ('Data_1/p0_1', 'p0')
save ('Data_1/parameters_1', 'par')

%% 3) Plot results \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -----------------------------------------------------

% Plot: Thigh angle estimate - acceleration-based,
%       KF, and EKF.
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

if tikz   
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
end
       
n = n + 1;

% Plot: Shank angle estimate - acceleration-based,
%       KF, and EKF.
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
    
if tikz   
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
end
        
n = n + 1;        

% Compute root-mean-square error.
% -Thigh
RMSE_acc = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
                - pitch_acc_right_thigh(n1:n2)').^2));       
RMSE_KF = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
               - pitch_KF_right_thigh(n1:n2)').^2));
RMSE_EKF = sqrt(mean((pitch_QS_right_thigh(n1:n2) ...
             - 90 - pitch_EKF_right_thigh(n1:n2)).^2));
RMSE = [RMSE_acc, RMSE_KF, RMSE_EKF];
% -Shank
RMSE_acc = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
                - pitch_acc_right_shank(n1:n2)').^2));       
RMSE_KF = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
               - pitch_KF_right_shank(n1:n2)').^2));
RMSE_EKF = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
             - 90 - pitch_EKF_right_thigh(n1:n2) ...
             - pitch_EKF_right_shank(n1:n2)).^2));
RMSE = [RMSE; RMSE_acc, RMSE_KF, RMSE_EKF;];

% Compute sum of seperate RMSEs.
RMSE = [RMSE; RMSE(1, 1) + RMSE(2, 1), RMSE(1, 2) + ...
        RMSE(2, 2), RMSE(1, 3) + RMSE(2, 3)];

% Plot bar graph and values on top of the bars.    
figure(n);
b = bar(RMSE, 0.3);
offset = 0.8;
yb = cat(1, b.YData); 
xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
hold on;

for i = 1:3  
   for j = 1:3
        text(xb(j, i), yb(j, i) + offset, ...
            ['\scriptsize ', num2str(RMSE(i, j), ...
            '$%0.2f$')], 'rotation', 0, ...
            'interpreter', 'latex', ...
            'HorizontalAlignment','center');
   end
end

b(1).FaceColor = [0.8500    0.3250    0.0980];
b(2).FaceColor = [0.9290    0.6940    0.1250];
b(3).FaceColor = [0.4940    0.1840    0.5560];

text(1.4, 18, ['\scriptsize $\frac{\operatorname', ...
        '{RMSE}_{EKF_{thigh}} + \operatorname', ...
        '{RMSE}_{EKF_{shank}}}{\operatorname', ...
        '{RMSE}_{KF_{thigh}} + \operatorname', ...
        '{RMSE}_{KF_{shank}}} = ', ...
        num2str((RMSE(1, 3) + RMSE(2, 3)) / ...
        (RMSE(1, 2) + RMSE(2, 2)), '%0.2f$')], ...
        'interpreter','latex');

ylim([0, max(max(RMSE)) + 2]);
ylabel('Root-mean-square error in $^{\circ}$', ...
       'interpreter','latex');

labels = {'Thigh', 'Shank', 'Thigh + Shank'};
format_ticks(gca, labels, [], [], [], 0);

legend('Acceleration-based', 'Kalman filter', ...
          'Extended Kalman filter');

if tikz   
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
end
      
n = n + 1;
         
% Plot: Acceleration-based pitch angle shank - corrected.
n1 = 4 * f + 1;
n2 = 20 * f;
figure(n);
hold on;
plot(time(n1:n2), pitch_QS_right_shank(n1:n2) - 90, ...
     'linewidth', 1);
plot(time(n1:n2), pitch_acc_right_shank(n1:n2) - 90);
plot(time(n1:n2), theta12_c(n1:n2), 'linewidth', 1);
 
xlabel('Time $t$ in s', 'interpreter','latex');
ylabel(['Pitch angle $\theta_1 + \theta_2$ in ', ...
        '$^{\circ}$'], 'interpreter', 'latex');
legend('Reference', 'Accelerometer-based', ...
       'Accelerometer based - corrected');
  
% Compute root-mean-square error.       
RMSE_acc = sqrt(mean((pitch_QS_right_shank(n1:n2) ...
                - pitch_acc_right_shank(n1:n2)').^2));       
RMSE_acc_corr = sqrt(mean((pitch_QS_right_shank(n1:n2)...
                - 90 - theta12_c(n1:n2)).^2));
RMSE_acc = [RMSE_acc, RMSE_acc_corr];
   
text(4.5,-130, ['\scriptsize $\operatorname{RMSE}', ...
                '_{Acc} = ', num2str(RMSE_acc(1), ...
                '%0.2f$')], 'interpreter','latex');
text(4.5,-140, ['\scriptsize $\operatorname{RMSE}', ...
                '_{Acc-corrected} = ', ...
                num2str(RMSE_acc(2), '%0.2f$')], ...
                'interpreter','latex');
text(4.5,-150, ['\scriptsize $\frac{\operatorname', ...
                '{RMSE}_{Acc-corrected}}', ...
                '{\operatorname{RMSE}_{Acc}} = ', ...
                num2str(RMSE_acc(2) / RMSE_acc(1), ...
                '%0.2f$')], 'interpreter','latex');
   
if tikz   
matlab2tikz(['../tikz/experiment_', num2str(n), ...
             '.tikz'], 'height', '\figureheight', ...
             'width', '\figurewidth');
end
%%         
% Display the results.
fprintf('-----------------Results---------------\n');
fprintf(['AB:\n RMSE_thigh: %0.4f\n RMSE_shank: ', ...
         '%0.4f\n RMSE_sum:   %0.4f\n\n'], RMSE(:, 1));
fprintf(['KF:\n RMSE_thigh: %0.4f\n RMSE_shank: ', ...
         '%0.4f\n RMSE_sum:   %0.4f\n\n'], RMSE(:, 2));
fprintf(['EKF:\n RMSE_thigh: %0.4f\n RMSE_shank: ', ...
         '%0.4f\n RMSE_sum:   %0.4f\n\n'], RMSE(:, 3));
fprintf(['Improvement:\n', ...
         'RMSE_SUM_EKF / RMSE_SUM_KF:  %0.4f\n\n'], ...
         RMSE(3, 3) / RMSE(3, 2));
fprintf(['Acceleration correction:\n', ...
         'RMSE_Acc / RMSE_Acc_corr:  %0.4f\n\n'], ...
         RMSE_acc(2) / RMSE_acc(1));
fprintf(['EKF parameters: \n -sigma_d_sq:  ', ...
         '%0.4f\n -sigma_t1_sq: ', ...
         '%0.4f\n -sigma_t2_sq: ', ...
         '%0.4f\n -sigma_b_sq:  ', ...
         '%0.4f\n -sigma_1_sq:  ', ...
         '%0.4f\n -sigma_2_sq:  ', ...
         '%0.4f\n -sigma_s3_sq: ', ...
         '%0.4f\n -sigma_f3_sq: ', ...
         '%0.4f\n -sigma_s4_sq: ', ...
         '%0.4f\n -sigma_f4_sq: ', ...
         '%0.4f\n'], par);
fprintf('---------------------------------------\n');

    