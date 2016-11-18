%% 1) Initialisation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% --------------------------------------------------------------------------------------------

clear all; close all; clc;

% Choose data sets for evaluation.
first_data_set = 1;
last_data_set = 40;

% Set length of thigh and shank.
l1 = 0.35;  % 0.35;  
l2 = 0.25;  % 0.25;

% Control plotting.
plot_raw_signals = false;             % Plot raw acceleration-based angle estimates before quad shift correction.
plot_KF_signals = false;              % Plot Kalman filter signals.
plot_acc_KF_signals = false;          % Plot Kalman filter signals including acceleration-based angle estimates.
plot_EKF_signals = false;             % Plot extended Kalman filter signals.
plot_KF_EKF_signals = true;           % Plot Kalman filter signals and extended Kalman filter signals.
plot_corracc_EKF_signals = false;     % Plot extended Kalman filter signals including motion-corrected acceleration-based angle estimates.
plot_acc_corracc_EKF_signals = true;  % Plot extended Kalman filter signals including acceleration-based and motion-corrected acceleration-based angle estimates.
plot_error_graphs = true;             % Plot root-mean-square errors.
plot_parameters_KF = false;            % Plot result of the parameter optimization

print_results = true;                 % Print individual results to the console.
save_figures = true;                 % Save figures to file.

% Control optimisation: if true no optimisation is carried out.
optimise_KF = false;
optimise_EKF = false;

% Save console output to file.
diary('console_output.txt');

% Set initial value of Kalman filter parameters;
measurement_noise_variance = 1.764798761889; %1000;
process_noise_variance = 0.000000319814; %0.001;

% Set initial values of parameters for optimisation of the extended Kalman filter.
sigma_s3 = 10;       % 3.5       measurement noise variance theta_1 - slow motion
sigma_s4 = 10;       % 30        measurement noise variance theta_1 + theta_2 - slow motion
sigma_f3 = 50;       % 30        measurement noise variance theta_1 - fast motion
sigma_f4 = 50;       % 300       measurement noise variance theta_1 + theta_2 - fast motion
sigma_b = 0.00001;   % 0.001     process noise variance gyroscope bias theta_1
sigma_t1 = 0.01;     % 0.05      process noise variance theta_1, maps to sigma_t2 in the internal of the EKF.
% sigma_d = 10;                 process noise variance of position along x and z-axis, defined internally in the EKF.

% Stores the initial params and the overall RMSE for optimisation.
load('params_results.mat')
params_results = [params_results; ...
    measurement_noise_variance, 0, 0, 0, 0, process_noise_variance, 0, 0, 0, 0, 1, -1, ...
    sigma_s3, sigma_s4, sigma_f3, sigma_f4, sigma_b, sigma_t1, 1, 100, 100];

%Initialise
%params_results = zeros(1, 21);
%save('params_results', 'params_results');

% Specify period in seconds of still standing at the beginning of the
% trial, which is used for the bias correction.
bias_corr_start = 1;
bias_corr_end = 2;

% Import GaitWatch and WaGyroMag functions library.
% All existing functions have to be called using
% either 'gw.functionName' or 'wag.functionName'.
gw = gwLibrary;
wag = wagLibrary;

% Set value of the magnitude of the gravity vector.
g = 9.81;

% Store file names in a cell array.
file_names = {'subj01_GaitWatch_1320_2 km.mat' ...
    'subj01_GaitWatch_1321_4 km.mat' ...
    'subj01_GaitWatch_1322_6 km.mat' ...
    'subj02_GaitWatch_1323_2 km.mat' ...
    'subj02_GaitWatch_1324_4 km.mat' ...
    'subj02_GaitWatch_1325_6 km.mat' ...
    'subj03_GaitWatch_1326_2 km.mat' ...
    'subj03_GaitWatch_1327_4 km.mat' ...
    'subj03_GaitWatch_1328_6 km.mat' ...
    'subj04_GaitWatch_1329_2 km.mat' ...
    'subj04_GaitWatch_1333_4 km.mat' ...
    ...  'subj04_GaitWatch_1334_6 km.mat' ...    Left shank, right thigh, right shank - strange quad shifts
    'subj04_GaitWatch_1336_6 km.mat' ...
    'subj05_GaitWatch_1337_2 km.mat' ...
    'subj05_GaitWatch_1338_2 km.mat' ...
    'subj05_GaitWatch_1339_4 km.mat' ...
    'subj05_GaitWatch_1340_6 km.mat' ...
    'subj06_GaitWatch_1341_2 km.mat' ...
    'subj06_GaitWatch_1342_4 km.mat' ...
    'subj06_GaitWatch_1343_6 km.mat' ...
    'subj07_GaitWatch_1344_2 km.mat' ...
    'subj07_GaitWatch_1345_2 km.mat' ...
    'subj07_GaitWatch_1346_4 km.mat' ...
    'subj07_GaitWatch_1347_6 km.mat' ...
    'subj08_GaitWatch_1348_2 km.mat' ...
    ...  'subj08_GaitWatch_1349_4 km.mat' ...    Left shank, left thigh - strange quad shifts
    ...  'subj08_GaitWatch_1350_6 km.mat' ...    Left shank, right shank - strange quad shifts
    'subj11_GaitWatch_1365_2 km.mat' ...
    'subj11_GaitWatch_1366_4 km.mat' ...
    ...  'subj11_GaitWatch_1367_6 km.mat' ...    Right thigh - strange quad shifts
    'subj12_GaitWatch_1358_2 km.mat' ...
    'subj12_GaitWatch_1359_4 km.mat' ...
    ...  'subj12_GaitWatch_1360_6 km.mat' ...    Left shank, right thigh, right shank - strange quad shifts
    'subj13_GaitWatch_1361_2 km.mat' ...
    'subj13_GaitWatch_1362_4 km.mat' ...
    'subj13_GaitWatch_1363_6 km.mat' ...
    'subj13_GaitWatch_1364_6 km.mat'};

%% 2) Evaluation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% --------------------------------------------------------------------------------------------

% Store the number of data sets.
n_data_sets = length(file_names);

% Initialise error matrix. It contains the errors of the thigh, shank, and
% thigh + shank angles determined acceleration-based, with the Kalman
% filter, and with the extended Kalman filter, respectively, in the order
% depicted in the error graphs.
error_matrix = zeros(n_data_sets, 28);

% Initialise parameter matrix. This will contain the 4 x 2 parameters alpha and
% beta of the Kalman filter as well as the 2 x 10 parameters of the extended Kalman
% filter, that is for the left and right body side, respectively.
parameter_matrix = zeros(n_data_sets, 28);

% Initialise vector that stores the walking speeds.
walking_speeds = zeros(n_data_sets, 1);

% Correct value of last_data_set if out of range.
if last_data_set > n_data_sets
    
    last_data_set = n_data_sets;
    
end

for index = first_data_set:last_data_set %n_data_sets
    
    file_name = file_names{index};
    
    % Load signals of the GaitWatch and the Qualisys motion capture system.
    load(file_name);
    
    % Store length of the current data set.
    len = length(g_Y_right_shank_1_C);
    
    % Extract subject number.
    subject = file_name(5:6);
    
    % Extract and store walking speed found in the file name.
    walking_speed = str2double(file_name(23));
    walking_speeds(index) = walking_speed;
    
    % Compute pitch using acceleration.
    pitch_acc_left_thigh = atan2d(a_Z_left_thigh_1_C, a_X_left_thigh_1_C);
    pitch_acc_left_shank = atan2d(a_Z_left_shank_1_C, a_X_left_shank_1_C);
    pitch_acc_right_thigh = atan2d(a_Z_right_thigh_1_C, a_X_right_thigh_1_C);
    pitch_acc_right_shank = atan2d(a_Z_right_shank_1_C, a_X_right_shank_1_C);
    
    % Adopt coordinate system where the horizontal equals 0 degrees.
    pitch_acc_left_thigh = pitch_acc_left_thigh - 180;
    pitch_acc_left_shank = pitch_acc_left_shank - 180;
    pitch_acc_right_thigh = pitch_acc_right_thigh - 180;
    pitch_acc_right_shank = pitch_acc_right_shank - 180;
    
    % Correct quadrant shifts.
    pitch_acc_left_thigh = gw.correct_quad_shifts(pitch_acc_left_thigh, 'deg');
    pitch_acc_left_shank = gw.correct_quad_shifts(pitch_acc_left_shank, 'deg');
    pitch_acc_right_thigh = gw.correct_quad_shifts(pitch_acc_right_thigh, 'deg');
    pitch_acc_right_shank = gw.correct_quad_shifts(pitch_acc_right_shank, 'deg');
    
    % Convert Rad to Deg and adopt coordinate system where the horizontal equals 0 degrees.
    pitch_QS_left_thigh = pitch_QS_left_thigh * 360 / (2 * pi) - 90;
    pitch_QS_left_shank = pitch_QS_left_shank * 360 / (2 * pi) - 90;
    pitch_QS_right_thigh = pitch_QS_right_thigh * 360 / (2 * pi) - 90;
    pitch_QS_right_shank = pitch_QS_right_shank * 360 / (2 * pi) - 90;
    
    % Multiply end of the bias correction period by sample frequency to get
    % an index instead of a time value.
    bias_corr_end = 2 * f;
    
    % Compute Qualisys biases.
    bias_left_thigh = mean(pitch_acc_left_thigh(bias_corr_start:bias_corr_end)) - mean(pitch_QS_left_thigh(bias_corr_start:bias_corr_end));
    bias_left_shank = mean(pitch_acc_left_shank(bias_corr_start:bias_corr_end)) - mean(pitch_QS_left_shank(bias_corr_start:bias_corr_end));
    bias_right_thigh = mean(pitch_acc_right_thigh(bias_corr_start:bias_corr_end)) - mean(pitch_QS_right_thigh(bias_corr_start:bias_corr_end));
    bias_right_shank = mean(pitch_acc_right_shank(bias_corr_start:bias_corr_end)) - mean(pitch_QS_right_shank(bias_corr_start:bias_corr_end));
    
    % Add bias.
    pitch_QS_left_thigh = pitch_QS_left_thigh + bias_left_thigh;
    pitch_QS_left_shank = pitch_QS_left_shank + bias_left_shank;
    pitch_QS_right_thigh = pitch_QS_right_thigh + bias_right_thigh;
    pitch_QS_right_shank = pitch_QS_right_shank + bias_right_shank;
    
    % Create time vector.
    time = (0:len - 1)' * 1 / f;
    
    
    %% 3) Optimisation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ----------------------------------------------------------------------------------------
    
    % Constitute initial paramter vector for optimisation of KF.
    p0_KF_left_thigh = [measurement_noise_variance process_noise_variance];
    p0_KF_left_shank = [measurement_noise_variance process_noise_variance];
    p0_KF_right_thigh = [measurement_noise_variance process_noise_variance];
    p0_KF_right_shank = [measurement_noise_variance process_noise_variance];
    
    
    if optimise_KF
        
        % Optimise filter parameters left thigh.
        
        % Call the optimization routine.
        [xmin, fmin, ct] = gw.optimize_KF(g_Y_left_thigh_1_C', pitch_acc_left_thigh', ...
            f, var(pitch_acc_left_thigh), var(pitch_acc_left_thigh), var(g_Y_left_thigh_1_C), ...
            pitch_acc_left_thigh(1), pitch_QS_left_thigh', p0_KF_left_thigh, 1);
        
        fprintf('\n------ KF OPTIMISATION left thigh ------\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf('Optimal parameters are: \n -Alpha: %0.4f \n -Beta: %0.4f\n', xmin(1), xmin(2));
        fprintf('----------------------------------------\n\n');
        
        % Extract optimal parameters.
        opt_alpha_KF_left_thigh = xmin(1);
        opt_beta_KF_left_thigh = xmin(2);
        
        % Optimise filter parameters left shank.
        
        % Call the optimization routine.
        [xmin, fmin, ct] = gw.optimize_KF(g_Y_left_shank_1_C', pitch_acc_left_shank', ...
            f, var(pitch_acc_left_shank), var(pitch_acc_left_shank), var(g_Y_left_shank_1_C), ...
            pitch_acc_left_shank(1), pitch_QS_left_shank', p0_KF_left_shank, 1);
        
        fprintf('\n------ KF OPTIMISATION left shank ------\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf('Optimal parameters are: \n -Alpha: %0.4f \n -Beta: %0.4f\n', xmin(1), xmin(2));
        fprintf('----------------------------------------\n\n');
        
        % Extract optimal parameters.
        opt_alpha_KF_left_shank = xmin(1);
        opt_beta_KF_left_shank = xmin(2);
        
        % Optimise filter parameters right thigh.
        
        % Call the optimization routine.
        [xmin, fmin, ct] = gw.optimize_KF(g_Y_right_thigh_1_C', pitch_acc_right_thigh', ...
            f, var(pitch_acc_right_thigh), var(pitch_acc_right_thigh), var(g_Y_right_thigh_1_C), ...
            pitch_acc_right_thigh(1), pitch_QS_right_thigh', p0_KF_right_thigh, 1);
        
        fprintf('\n------ KF OPTIMISATION right thigh ------\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf('Optimal parameters are: \n -Alpha: %0.4f \n -Beta: %0.4f\n', xmin(1), xmin(2));
        fprintf('----------------------------------------\n\n');
        
        % Extract optimal parameters.
        opt_alpha_KF_right_thigh = xmin(1);
        opt_beta_KF_right_thigh = xmin(2);
        
        % Optimise filter parameters right shank.
        
        % Call the optimization routine.
        [xmin, fmin, ct] = gw.optimize_KF(g_Y_right_shank_1_C', pitch_acc_right_shank', ...
            f, var(pitch_acc_right_shank), var(pitch_acc_right_shank), var(g_Y_right_shank_1_C), ...
            pitch_acc_right_shank(1), pitch_QS_right_shank', p0_KF_right_shank, 1);
        
        fprintf('\n------ KF OPTIMISATION right shank ------\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf('Optimal parameters are: \n -Alpha: %0.4f \n -Beta: %0.4f\n', xmin(1), xmin(2));
        fprintf('----------------------------------------\n\n');
        
        % Extract optimal parameters.
        opt_alpha_KF_right_shank = xmin(1);
        opt_beta_KF_right_shank = xmin(2);
        
    else
        
        % If no optimisation is carried out use initial parameters
        opt_alpha_KF_left_thigh = p0_KF_left_thigh(1);
        opt_beta_KF_left_thigh = p0_KF_left_thigh(2);
        
        opt_alpha_KF_left_shank = p0_KF_left_shank(1);
        opt_beta_KF_left_shank = p0_KF_left_shank(2);
        
        opt_alpha_KF_right_thigh = p0_KF_right_thigh(1);
        opt_beta_KF_right_thigh = p0_KF_right_thigh(2);
        
        opt_alpha_KF_right_shank = p0_KF_right_shank(1);
        opt_beta_KF_right_shank = p0_KF_right_shank(2);
        
    end
    
    % Constitute initial paramter vector for optimisation of EKF.
    p0_EKF_left = params_results(end, 13:18);
    p0_EKF_right = params_results(end, 13:18);
    
    if optimise_EKF
        
        % Optimise left EKF.
        
        % Call the optimization routine.
        [p_opt_left, fmin, ct] = gw.optimize_EKF(g_Y_left_thigh_1_C', g_Y_left_shank_1_C', ...
            a_X_left_thigh_1_C', a_Z_left_thigh_1_C', a_X_left_shank_1_C', a_Z_left_shank_1_C', ...
            f, l1, l2, [pitch_QS_left_thigh; pitch_QS_left_shank], p0_EKF_left, 1);
        
        fprintf('\n-- EKF OPTIMISATION left thigh + shank -\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf(['Optimal parameters are: \n -sigma_s3: ', ...
            '%0.4f\n -sigma_s4: %0.4f\n -sigma_f3: ', ...
            '%0.4f\n -sigma_f4: %0.4f\n -sigma_b: ', ...
            '%0.4f\n -sigma_t1: %0.4f'], p_opt_left);
        fprintf('\n----------------------------------------\n\n')
        
        % Optimise right EKF.
        
        % Call the optimization routine.
        [p_opt_right, fmin, ct] = gw.optimize_EKF(g_Y_right_thigh_1_C', g_Y_right_shank_1_C', ...
            a_X_right_thigh_1_C', a_Z_right_thigh_1_C', a_X_right_shank_1_C', a_Z_right_shank_1_C', ...
            f, l1, l2, [pitch_QS_right_thigh; pitch_QS_right_shank], p0_EKF_right, 1);
        
        fprintf('\n-- EKF OPTIMISATION right thigh + shank -\n');
        fprintf('The optimisation process finished in %d iterations.\n', ct);
        fprintf('The minimum RMSE found is: %0.4f\n', fmin);
        fprintf(['Optimal parameters are: \n -sigma_s3: ', ...
            '%0.4f\n -sigma_s4: %0.4f\n -sigma_f3: ', ...
            '%0.4f\n -sigma_f4: %0.4f\n -sigma_b: ', ...
            '%0.4f\n -sigma_t1: %0.4f'], p_opt_right);
        fprintf('\n----------------------------------------\n\n')
        
    else
        
        p_opt_left = p0_EKF_left;
        p_opt_right = p0_EKF_right;
        
    end
    
    %% 4) Apply Kalman filters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ----------------------------------------------------------------------------------------
    
    % Kalman filter
    pitch_KF_left_thigh = gw.fusion_KF(g_Y_left_thigh_1_C, pitch_acc_left_thigh, ...
        f, var(pitch_acc_left_thigh), var(pitch_acc_left_thigh), var(g_Y_left_thigh_1_C), ...
        opt_alpha_KF_left_thigh, opt_beta_KF_left_thigh, pitch_acc_left_thigh(1));
    
    pitch_KF_left_shank = gw.fusion_KF(g_Y_left_shank_1_C, pitch_acc_left_shank, ...
        f, var(pitch_acc_left_shank), var(pitch_acc_left_shank), var(g_Y_left_shank_1_C), ...
        opt_alpha_KF_left_shank, opt_beta_KF_left_shank, pitch_acc_left_shank(1));
    
    pitch_KF_right_thigh = gw.fusion_KF(g_Y_right_thigh_1_C, pitch_acc_right_thigh, ...
        f, var(pitch_acc_right_thigh), var(pitch_acc_right_thigh), var(g_Y_right_thigh_1_C), ...
        opt_alpha_KF_right_thigh, opt_beta_KF_right_thigh, pitch_acc_right_thigh(1));
    
    pitch_KF_right_shank = gw.fusion_KF(g_Y_right_shank_1_C, pitch_acc_right_shank, ...
        f, var(pitch_acc_right_shank), var(pitch_acc_right_shank), var(g_Y_right_shank_1_C), ...
        opt_alpha_KF_right_shank, opt_beta_KF_right_shank, pitch_acc_right_shank(1));
    
    % Store parameters to parameter matrix.
    parameter_matrix(index, 1:8) = [opt_alpha_KF_left_thigh, opt_beta_KF_left_thigh, ...
        opt_alpha_KF_left_shank, opt_beta_KF_left_shank, ...
        opt_alpha_KF_right_thigh, opt_beta_KF_right_thigh, ...
        opt_alpha_KF_right_shank, opt_beta_KF_right_shank];
    
    params_results(end, 2:5) = [opt_alpha_KF_left_thigh, opt_alpha_KF_left_shank, opt_alpha_KF_right_thigh, opt_alpha_KF_right_shank];
    params_results(end, 7:10) = [opt_beta_KF_left_thigh, opt_beta_KF_left_shank, opt_beta_KF_right_thigh, opt_beta_KF_right_shank];
    
    
    % Extended Kalman filter.
    % Additionally, store the internal state vector at each time step in x, the motion based
    % acceleration in a_m, and the angle estimate theta_1 + theta_2, based on the corrected
    % acceleration signal in theta12_c.
    [pitch_EKF_left_thigh, pitch_EKF_left_shank, pitch_acc_c_left_shank, a_m_left, x_left, par_left] ...
        = fusion_EKF(g_Y_left_thigh_1_C', g_Y_left_shank_1_C', a_X_left_thigh_1_C', ...
        a_Z_left_thigh_1_C', a_X_left_shank_1_C', a_Z_left_shank_1_C', f, l1, l2, p_opt_left);
    
    [pitch_EKF_right_thigh, pitch_EKF_right_shank, pitch_acc_c_right_shank, a_m_right, x_right, par_right] ...
        = fusion_EKF(g_Y_right_thigh_1_C', g_Y_right_shank_1_C', a_X_right_thigh_1_C', ...
        a_Z_right_thigh_1_C', a_X_right_shank_1_C', a_Z_right_shank_1_C', f, l1, l2, p_opt_right);
    
    % Transpose and consider that the shank angle in the EKF is measured
    % with respect to the thigh. That is it needs to be added in order to
    % obtain the shank angle with respect to the horizontal.
    pitch_EKF_left_thigh = pitch_EKF_left_thigh';
    pitch_EKF_left_shank = pitch_EKF_left_thigh + pitch_EKF_left_shank';
    
    pitch_EKF_right_thigh = pitch_EKF_right_thigh';
    pitch_EKF_right_shank = pitch_EKF_right_thigh + pitch_EKF_right_shank';
    
    % Store parameters in parameter matrix.
    parameter_matrix(index, 9:28) = [par_left, par_right];
    
    %% 5) Compute Errors \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ----------------------------------------------------------------------------------------
    
    % Compute root-mean-square errors of Kalman filters left.
    RMSE_acc_left_thigh = sqrt(mean((pitch_QS_left_thigh - pitch_acc_left_thigh).^2));
    RMSE_KF_left_thigh = sqrt(mean((pitch_QS_left_thigh - pitch_KF_left_thigh).^2));
    RMSE_EKF_left_thigh = sqrt(mean((pitch_QS_left_thigh - pitch_EKF_left_thigh).^2));
    RMSE_left_thigh = [RMSE_acc_left_thigh, RMSE_KF_left_thigh, RMSE_EKF_left_thigh];
    
    RMSE_acc_left_shank = sqrt(mean((pitch_QS_left_shank - pitch_acc_left_shank).^2));
    RMSE_KF_left_shank = sqrt(mean((pitch_QS_left_shank - pitch_KF_left_shank).^2));
    RMSE_EKF_left_shank = sqrt(mean((pitch_QS_left_shank - pitch_EKF_left_shank).^2));
    RMSE_left_shank = [RMSE_acc_left_shank, RMSE_KF_left_shank, RMSE_EKF_left_shank];
    
    RMSE_left_sum = RMSE_left_thigh + RMSE_left_shank;
    
    % Compute root-mean-square error of motion-based acceleration correction.
    RMSE_acc_corr_left_shank = sqrt(mean((pitch_QS_left_shank - pitch_acc_c_left_shank').^2));
    
    RMSE_left = [RMSE_left_thigh; RMSE_left_shank; RMSE_left_sum];
    
    % Compute root-mean-square errors of Kalman filters right.
    RMSE_acc_right_thigh = sqrt(mean((pitch_QS_right_thigh - pitch_acc_right_thigh).^2));
    RMSE_KF_right_thigh = sqrt(mean((pitch_QS_right_thigh - pitch_KF_right_thigh).^2));
    RMSE_EKF_right_thigh = sqrt(mean((pitch_QS_right_thigh - pitch_EKF_right_thigh).^2));
    RMSE_right_thigh = [RMSE_acc_right_thigh, RMSE_KF_right_thigh, RMSE_EKF_right_thigh];
    
    RMSE_acc_right_shank = sqrt(mean((pitch_QS_right_shank - pitch_acc_right_shank).^2));
    RMSE_KF_right_shank = sqrt(mean((pitch_QS_right_shank - pitch_KF_right_shank).^2));
    RMSE_EKF_right_shank = sqrt(mean((pitch_QS_right_shank - pitch_EKF_right_shank).^2));
    RMSE_right_shank = [RMSE_acc_right_shank, RMSE_KF_right_shank, RMSE_EKF_right_shank];
    
    RMSE_right_sum = RMSE_right_thigh + RMSE_right_shank;
    
    % Compute root-mean-square error of motion-based acceleration correction.
    RMSE_acc_corr_right_shank = sqrt(mean((pitch_QS_right_shank - pitch_acc_c_right_shank').^2));
    
    RMSE_right = [RMSE_right_thigh; RMSE_right_shank; RMSE_right_sum];
    
    % Compute sum of RMSE of the shank, raw acceleration based and after
    % acceleration correction.
    RMSE_acc_shank_sum = RMSE_acc_left_shank + RMSE_acc_right_shank;
    RMSE_acc_corr_shank_sum = RMSE_acc_corr_left_shank + RMSE_acc_corr_right_shank;
    
    
    
    % Store all RMSE in one matrix.
    RMSE = [RMSE_left; RMSE_right; RMSE_left + RMSE_right; RMSE_acc_shank_sum, RMSE_acc_corr_shank_sum, 0];
    
    error_matrix(index, :) = [RMSE_left_thigh, RMSE_left_shank, RMSE_left_sum, ...
        RMSE_right_thigh, RMSE_right_shank, RMSE_right_sum, ...
        RMSE_left_thigh + RMSE_right_thigh, RMSE_left_shank + RMSE_right_shank, ...
        RMSE_left_thigh + RMSE_right_thigh + RMSE_left_shank + RMSE_right_shank, ...
        RMSE_acc_corr_shank_sum];
    
    %% 6) Plots \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ----------------------------------------------------------------------------------------
    
    if plot_raw_signals
        
        % Plot raw signals without quad shift correction.
        
        % Left thigh
        figure();
        hold on;
        plot(time, pitch_QS_left_thigh, 'linewidth', 2);
        plot(time, pitch_acc_left_thigh, 'linewidth', 0.5);
        plot(0);
        plot(0);
        title(['Left thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(time, pitch_acc_left_shank, 'linewidth', 0.5);
        plot(0);
        plot(0);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right thigh
        figure();
        hold on;
        plot(time, pitch_QS_right_thigh, 'linewidth', 2);
        plot(time, pitch_acc_right_thigh, 'linewidth', 0.5);
        plot(0);
        plot(0);
        title(['Right thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(time, pitch_acc_right_shank, 'linewidth', 0.5);
        plot(0);
        plot(0);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_KF_signals
        
        % Plot Kalman filter signals.
        
        % Left thigh
        figure();
        hold on;
        plot(time, pitch_QS_left_thigh, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_left_thigh, 'linewidth', 1);
        plot(0);
        title(['Left thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_left_shank, 'linewidth', 1);
        plot(0);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right thigh
        figure();
        hold on;
        plot(time, pitch_QS_right_thigh, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_right_thigh, 'linewidth', 1);
        plot(0);
        title(['Right thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_right_shank, 'linewidth', 1);
        plot(0);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_acc_KF_signals
        
        % Plot Kalman filter signals.
        
        % Left thigh
        figure();
        hold on;
        plot(time, pitch_QS_left_thigh, 'linewidth', 2);
        plot(time, pitch_acc_left_thigh, 'linewidth', 0.5);
        plot(time, pitch_KF_left_thigh, 'linewidth', 1);
        plot(0);
        title(['Left thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(time, pitch_acc_left_shank, 'linewidth', 0.5);
        plot(time, pitch_KF_left_shank, 'linewidth', 1);
        plot(0);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right thigh
        figure();
        hold on;
        plot(time, pitch_QS_right_thigh, 'linewidth', 2);
        plot(time, pitch_acc_right_thigh, 'linewidth', 0.5);
        plot(time, pitch_KF_right_thigh, 'linewidth', 1);
        plot(0);
        title(['Right thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(time, pitch_acc_right_shank, 'linewidth', 0.5);
        plot(time, pitch_KF_right_shank, 'linewidth', 1);
        plot(0);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_EKF_signals
        
        % Plot extended Kalman filter signals.
        
        % Left thigh
        figure();
        hold on;
        plot(time, pitch_QS_left_thigh, 'linewidth', 2);
        plot(0);
        plot(0);
        plot(time, pitch_EKF_left_thigh, 'linewidth', 1);
        title(['Left thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(0);
        plot(0);
        plot(time, pitch_EKF_left_shank, 'linewidth', 1);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right thigh
        figure();
        hold on;
        plot(time, pitch_QS_right_thigh, 'linewidth', 2);
        plot(0);
        plot(0);
        plot(time, pitch_EKF_right_thigh, 'linewidth', 1);
        title(['Right thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(0);
        plot(0);
        plot(time, pitch_EKF_right_shank, 'linewidth', 1);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_acc_corracc_EKF_signals
        
        % Plot extended Kalman filter signals including acceleration-based and motion-corrected acceleration-based angle estimates.
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(time, pitch_acc_left_shank, 'linewidth', 0.5);
        plot(time, pitch_acc_c_left_shank, 'linewidth', 0.5);
        plot(time, pitch_EKF_left_shank, 'linewidth', 1);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Corrected acceleration-based', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(time, pitch_acc_right_shank, 'linewidth', 0.5);
        plot(time, pitch_acc_c_right_shank, 'linewidth', 0.5);
        plot(time, pitch_EKF_right_shank, 'linewidth', 1);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Corrected acceleration-based', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_corracc_EKF_signals
        
        % Plot extended Kalman filter signals including motion-corrected acceleration-based angle estimates.
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(time, pitch_acc_c_left_shank, 'linewidth', 0.5);
        plot(0);
        plot(time, pitch_EKF_left_shank, 'linewidth', 1);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Corrected acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(time, pitch_acc_c_right_shank, 'linewidth', 0.5);
        plot(0);
        plot(time, pitch_EKF_right_shank, 'linewidth', 1);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Corrected acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_KF_EKF_signals
        
        % Plot Kalman filter signals and extended Kalman filter signals.
        
        % Left thigh
        figure();
        hold on;
        plot(time, pitch_QS_left_thigh, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_left_thigh, 'linewidth', 1);
        plot(time, pitch_EKF_left_thigh, 'linewidth', 1);
        title(['Left thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Left shank
        figure();
        hold on;
        plot(time, pitch_QS_left_shank, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_left_shank, 'linewidth', 1);
        plot(time, pitch_EKF_left_shank, 'linewidth', 1);
        title(['Left shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right thigh
        figure();
        hold on;
        plot(time, pitch_QS_right_thigh, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_right_thigh, 'linewidth', 1);
        plot(time, pitch_EKF_right_thigh, 'linewidth', 1);
        title(['Right thigh - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
        % Right shank
        figure();
        hold on;
        plot(time, pitch_QS_right_shank, 'linewidth', 2);
        plot(0);
        plot(time, pitch_KF_right_shank, 'linewidth', 1);
        plot(time, pitch_EKF_right_shank, 'linewidth', 1);
        title(['Right shank - subject ', subject, ', walking speed ', num2str(walking_speed), ' km/h']);
        xlabel('Time $t$ in s');
        ylabel('Pitch angle $\theta_1$ in $^{\circ}$');
        legend('Qualisys', 'Acceleration-based', 'Kalman Filter', 'Extended Kalman Filter','Location', 'northwest');
        legend('boxoff');
        
    end
    
    if plot_error_graphs
        
        % Plot bar graph and values on top of the bars for left body side.
        figure();
        b = bar(RMSE, 0.3);
        % Offset for values on top of bars.
        offset = max(max(RMSE)) / 40;
        yb = cat(1, b.YData);
        xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
        hold on;
        
        % Add numbers on top of bars.
        for i = 1:9
            for j = 1:3
                text(xb(j, i), yb(j, i) + offset, [' ', num2str(RMSE(i, j), '$%0.2f$')], ...
                    'HorizontalAlignment','center');
            end
        end
        
        % Add numbers on top of bars of the last two bars.
        text(xb(1, 10), yb(1, 10) + offset, [' ', num2str(RMSE(10, 2), '$%0.2f$')], ...
                    'HorizontalAlignment','center');
        text(xb(2, 10), yb(2, 10) + offset, [' ', num2str(RMSE(10, 3), '$%0.2f$')], ...
                    'HorizontalAlignment','center');
        
        b(1).FaceColor = [0.8500    0.3250    0.0980];
        b(2).FaceColor = [0.9290    0.6940    0.1250];
        b(3).FaceColor = [0.4940    0.1840    0.5560];
        
        ylim([0, max(max(RMSE)) * 1.1 + offset]);
        ylabel('Root-mean-square error in $^{\circ}$');
        
        labels = {'Left THIGH', 'Left SHANK', 'Left THIGH + SHANK', ...
            'Right THIGH', 'Right SHANK', 'Right THIGH + SHANK', ...
            'Both THIGH', 'Both SHANK', 'Both THIGH + SHANK', 'Acceleration correction SHANK'};
        bx = gca;
        bx.XTickLabel = labels;
        bx.XTickLabelRotation = 45;
        
        title(['Errors - subject ', subject, ', walking speed ', ...
            num2str(walking_speed), ' km/h']);
        legend('Acceleration-based', 'Kalman filter', 'Extended Kalman filter', 'Location','northwest');
        legend('boxoff');
        
    end
    
    
    %% 7) Print results \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ----------------------------------------------------------------------------------------
    
    if print_results
        
        fprintf('-----------------------------------------------------------------------------\n');
        fprintf(['///// %i. Results subject ', subject, ', walking speed ', ...
            num2str(walking_speed), ' km/h /////\n'], index);
        fprintf('-----------------------------------------------------------------------------\n\n');
        
        fprintf('----------------- Results left body side ---------------\n');
        fprintf(['Acceleration-based:\n  RMSE_left_thigh: %0.4f\n  RMSE_left_shank: ', ...
            '%0.4f\n  RMSE_left_sum:   %0.4f\n\n'], RMSE_left(:, 1));
        fprintf(['Kalman filter:\n  RMSE_left_thigh: %0.4f\n  RMSE_left_shank: ', ...
            '%0.4f\n  RMSE_left_sum:   %0.4f\n\n'], RMSE_left(:, 2));
        fprintf(['Extended Kalman filter:\n  RMSE_left_thigh: %0.4f\n  RMSE_left_shank: ', ...
            '%0.4f\n  RMSE_left_sum:   %0.4f\n\n'], RMSE_left(:, 3));
        
        fprintf('Improvement:\n  1 - RMSE_left_EKF / RMSE_left_KF = %0.4f %%\n\n', ...
            (1 - RMSE_left_sum(3) / RMSE_left_sum(2)) * 100);
        
        fprintf('  1 - RMSE_acc_corr_left_shank / RMSE_acc_left_shank = %0.4f %%\n\n', ...
            (1 - RMSE_acc_corr_left_shank / RMSE_acc_left_shank) * 100);
        
        fprintf(['EKF parameters left : \n -sigma_d:  ', ...
            '%0.4f\n sigma_t1: ', ...
            '%0.4f\n sigma_t2: ', ...
            '%0.4f\n sigma_b:  ', ...
            '%0.4f\n sigma_1:  ', ...
            '%0.4f\n sigma_2:  ', ...
            '%0.4f\n sigma_s3: ', ...
            '%0.4f\n sigma_f3: ', ...
            '%0.4f\n sigma_s4: ', ...
            '%0.4f\n sigma_f4: ', ...
            '%0.4f\n\n\n'], par_left);
        
        fprintf('----------------- Results right body side ---------------\n');
        fprintf(['Acceleration-based:\n  RMSE_right_thigh: %0.4f\n  RMSE_right_shank: ', ...
            '%0.4f\n  RMSE_right_sum:   %0.4f\n\n'], RMSE_right(:, 1));
        fprintf(['Kalman filter:\n  RMSE_right_thigh: %0.4f\n  RMSE_right_shank: ', ...
            '%0.4f\n  RMSE_right_sum:   %0.4f\n\n'], RMSE_right(:, 2));
        fprintf(['Extended Kalman filter:\n  RMSE_right_thigh: %0.4f\n  RMSE_right_shank: ', ...
            '%0.4f\n  RMSE_right_sum:   %0.4f\n\n'], RMSE_right(:, 3));
        
        fprintf('Improvement:\n  1 - RMSE_right_EKF / RMSE_right_KF = %0.4f %%\n\n', ...
            (1 - RMSE_right_sum(3) / RMSE_right_sum(2)) * 100);
        
        fprintf('  1 - RMSE_acc_corr_right_shank / RMSE_acc_right_shank = %0.4f %%\n\n', ...
            (1 - RMSE_acc_corr_right_shank / RMSE_acc_right_shank) * 100);
        
        fprintf(['EKF parameters right : \n -sigma_d:  ', ...
            '%0.4f\n sigma_t1: ', ...
            '%0.4f\n sigma_t2: ', ...
            '%0.4f\n sigma_b:  ', ...
            '%0.4f\n sigma_1:  ', ...
            '%0.4f\n sigma_2:  ', ...
            '%0.4f\n sigma_s3: ', ...
            '%0.4f\n sigma_f3: ', ...
            '%0.4f\n sigma_s4: ', ...
            '%0.4f\n sigma_f4: ', ...
            '%0.4f\n\n\n'], par_right);
        
        fprintf('----------------- Results both body side ---------------\n');
        fprintf(['Acceleration-based:\n  RMSE_both_thigh: %0.4f\n  RMSE_both_shank: ', ...
            '%0.4f\n  RMSE_both_sum:   %0.4f\n\n'], RMSE_left(:, 1) + RMSE_right(:, 1));
        fprintf(['Kalman filter:\n  RMSE_both_thigh: %0.4f\n  RMSE_both_shank: ', ...
            '%0.4f\n  RMSE_both_sum:   %0.4f\n\n'], RMSE_left(:, 2) + RMSE_right(:, 2));
        fprintf(['Extended Kalman filter:\n  RMSE_both_thigh: %0.4f\n  RMSE_both_shank: ', ...
            '%0.4f\n  RMSE_both_sum:   %0.4f\n\n'], RMSE_left(:, 3) + RMSE_right(:, 3));
        
        fprintf('Improvement:\n  1 - RMSE_both_EKF / RMSE_both_KF = %0.4f %%\n\n', ...
            (1 - RMSE(9, 3) / RMSE(9, 2)) * 100);
        
        fprintf('  1 - RMSE_acc_corr_both_shank / RMSE_acc_both_shank = %0.4f %%\n\n', ...
            (1 - RMSE_acc_corr_shank_sum / RMSE_acc_shank_sum) * 100);
        
    end
    
end

% Combine results.
results = [sum(error_matrix(:, 22)), sum(error_matrix(:, 19)); ...
    sum(error_matrix(:, 28)), 0; ...
    sum(error_matrix(:, 23)), sum(error_matrix(:, 20)); ...
    sum(error_matrix(:, 24)), sum(error_matrix(:, 21)); ...
    ];

figure();
b = bar(results, 0.3, 'stacked'); 
% Offset for values on top of bars.
offset = max(max(results)) / 40;
yb = cat(1, b.YData);
xb = bsxfun(@plus, b(1).XData, [b.XOffset]');
hold on;

ylabel('Root-mean-square error in $^{\circ}$');
labels = {'Acceleration-based THIGH + SHANK', 'Acceleration-based corrected shank', 'Kalman filter THIGH + SHANK', 'Extended Kalman filter THIGH + SHANK'};
bx = gca;
bx.XTickLabel = labels;
bx.XTickLabelRotation = 0;
legend('Shank', 'Thigh');

title('Overall results');

if plot_parameters_KF
    
    data_set_no = 1:n_data_sets;
    
    % Scatter plot the filter parameters alpha and beta of the KF including the mean.
    figure();
    
    subplot(2,2,1);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 1), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 1), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 1), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 1));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 1));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 1));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Alpha left thigh');
    
    subplot(2,2,3);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 3), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 3), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 3), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 3));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 3));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 3));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Alpha left shank');
    
    subplot(2,2,2);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 5), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 5), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 5), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 5));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 5));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 5));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Alpha right thigh');
    
    subplot(2,2,4);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 7), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 7), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 7), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 7));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 7));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 7));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Alpha right shank');
    
    % Betas
    figure();
    
    subplot(2,2,1);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 2), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 2), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 2), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 2));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 2));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 2));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Beta left thigh');
    
    subplot(2,2,3);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 4), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 4), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 4), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 4));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 4));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 4));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Beta left shank');
    
    subplot(2,2,2);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 6), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 6), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 6), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 6));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 6));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 6));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Beta right thigh');
    
    subplot(2,2,4);
    hold on;
    scatter(data_set_no(walking_speeds == 2), parameter_matrix(walking_speeds == 2, 8), [], 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
    scatter(data_set_no(walking_speeds == 4), parameter_matrix(walking_speeds == 4, 8), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    scatter(data_set_no(walking_speeds == 6), parameter_matrix(walking_speeds == 6, 8), [], 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');
    legend('2 km/h', '4 km/h', '6 km/h');
    % Plot mean.
    mean_2 = mean(parameter_matrix(walking_speeds == 2, 8));
    mean_4 = mean(parameter_matrix(walking_speeds == 4, 8));
    mean_6 = mean(parameter_matrix(walking_speeds == 6, 8));
    hold on;
    plot([1 n_data_sets],[mean_2 mean_2], 'color', 'blue');
    plot([1 n_data_sets],[mean_4 mean_4], 'color', 'red');
    plot([1 n_data_sets],[mean_6 mean_6], 'color', 'green');
    title('Beta right shank');
    
    % Add mean of all parameters.
    
end

% Print end results
fprintf('\n\n\n-----------------------------------------------------------------------------\n');
fprintf('//////// OVERALL RESULTS ////////\n', index);
fprintf('-----------------------------------------------------------------------------\n\n');

fprintf('Acceleration-based:\n RMSE = %0.4f\n\n', sum(results(1, :)));
fprintf('Kalman filter:\n RMSE = %0.4f\n\n', sum(results(3, :)));
fprintf('Extended Kalman filter:\n RMSE = %0.4f\n\n', sum(results(4, :)));

fprintf('Improvement:\n  1 - RMSE_EKF / RMSE_KF = %0.4f %%\n\n', ...
    (1 - sum(results(4, :)) / sum(results(3, :))) * 100);

fprintf('  1 - RMSE_acc_corr_both_shank / RMSE_acc_both_shank = %0.4f %%\n\n', ...
    (1 - RMSE_acc_corr_shank_sum / RMSE_acc_shank_sum) * 100);


% Store end result of KF and EKF.
params_results(end, 11) = sum(results(3, :));
params_results(end, 19) = sum(results(4, :));

params_results(end, 20) = (1 - params_results(end, 19) / params_results(end, 11)) * 100;
params_results(end, 21) = (1 - RMSE_acc_corr_shank_sum / RMSE_acc_shank_sum) * 100;

% Save results.
save('results', 'parameter_matrix', 'error_matrix', 'results');
save('params_results', 'params_results');

% Save figures.
figHandles = get(0,'Children');
if save_figures
    savefig(figHandles, 'figures');
end


%%

%figure();
%stem(parameter_matrix(:, 23));

%figure();
%stem(error_matrix(:, end));
