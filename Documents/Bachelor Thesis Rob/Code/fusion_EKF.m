function [theta1, theta2, theta12_c, a_m, X] = ...
            fusion_EKF(gyro_thigh_y, gyro_shank_y, ...
                       acc_thigh_x, acc_thigh_z, ...
                       acc_shank_x, acc_shank_z, ...
                       fs, l1, l2)

% FUNCTION fusion_EKF applies an extended Kalman filter
% in order to fuse the accelerometer and gyroscope data
% and thus obtain an accurate orientation estimate of 
% the thighs and shanks.
%
% Input arguments:
% |_ 'gyro_thigh_y':  Row vector containing the angular 
%                     rate of the thigh about the 
%                     y-axis in radians per second.
% |_ 'gyro_shank_y':  Row vector containing the angular 
%                     rate of the shank about the 
%                     y-axis in radians per second.
% |_ 'acc_thigh_x':   Row vector containing the linear
%                     acceleration of the thigh along
%                     the x-axis in g.
% |_ 'acc_thigh_z':   Row vector containing the linear
%                     acceleration of the thigh along
%                     the z-axis in g.
% |_ 'acc_shank_x':   Row vector containing the linear
%                     acceleration of the shank along
%                     the x-axis in g.
% |_ 'acc_shank_z':   Row vector containing the linear
%                     acceleration of the shank along
%                     the z-axis in g.
% |_ 'fs':            Sampling frecuency in Hertz. Must
%                     be real positive.
% |_ 'l1':            Length of the thigh in m. Must be
%                     real positive.
% |_ 'l2':            Length of the shank in m. Must be
%                     real positive.
%
% Output:
% |_ 'theta1':        Row vector containing the thigh 
%                     angle with respect to the x-axis 
%                     of the world frame in radians.
% |_ 'theta2':        Row vector containing the shank 
%                     angle with respect to the thigh
%                     in radians.
% |_ 'theta2':        Row vector containing the shank 
%                     angle with respect to the thigh
%                     in radians.
% |_ 'a_m':           Row vector containing the 
%                     estimate of the acceleration that
%                     sensor 2 will see due to motion.
% |_ 'X':             10 x length(gyro_thigh_y) matrix
%                     containing the internal state
%                     vector at each instant as columns.
%
% IMPORTANT NOTE:     gyro_thigh_y, gyro_shank_y, 
%                     acc_thigh_x, acc_thigh_z, 
%                     acc_shank_x, and acc_shank_z 
%                     must have the same length. 
%                     Otherwise, an error will
%                     be returned.
% -----------------------------------------------------
% Authors:            Robin Weiss
% Entity:             University of Applied Sciences
%                     Munster, Munster, Germany
% Last modification:  13/05/2015
% -----------------------------------------------------

% 1) Check input arguments.
if ~isequal(length(gyro_thigh_y), ...
            length(gyro_shank_y), ...
            length(acc_thigh_x), ...
            length(acc_thigh_z), ...
            length(acc_shank_x), ...
            length(acc_shank_z))
    error(['Input arguments ''gyro_thigh_y'', ', ... 
           '''acc_thigh_x'', ''acc_thigh_z'', ', ...
           '''acc_shank_x'', ''acc_shank_z'', ', ...
           'must have the same length.']);
end

if (fs <= 0 || ~isreal(fs))
    error(['Input argument ''fs'' must be real', ...
           'positive.']);
end

if (l1 <= 0 || ~isreal(l1))
    error(['Input argument ''a1'' must be real', ...
          ' positive.']);
end

if (l2 <= 0 || ~isreal(l2))
    error(['Input parameter ''a2'' must be real', ...
           ' positive.']);
end

% 2) Import GaitWatch and WaGyroMag functions library.
%    All existing functions have to be called using
%    either 'gw.functionName' or 'wag.functionName'.
gw = gwLibrary;
wag = wagLibrary;

% 3) Compute the sampling period and the length of the
%    signal vectors.
Ts = 1 / fs;
len = length(gyro_thigh_y);

% 4) Set the magnitude of gravity.
gravity = 9.81;

% 5) Compute correction factor. This factor will be 
%    used throughout the entire code in order to
%    convert degrees into radians.
c_w = pi / 180;

% 6) Initialise output vectors.
theta1 = zeros(1, len);
theta2 = zeros(1, len);
theta12_c = zeros(1, len);
a_m = zeros(3, len);
X = zeros(10, len);

% 7) Compute intensity level.
lwin_ltsd = 20;    
threshold_ltsd = 4;    
shift_ltsd = 19;    
input_signal = sqrt(acc_shank_x.^2+acc_shank_z.^2);
[V_fsd, T_fsd] = wag.ltsd(input_signal', lwin_ltsd, ...
                      shift_ltsd, 512, threshold_ltsd);
                   
% 8) Determine marker signal.
[marker, ~] = gw.compEstMark(V_fsd, T_fsd, ...
                           input_signal, lwin_ltsd, ...
                           shift_ltsd);
                            
% INITIALISATION OF PARAMETERS %
                            
% 9) Compute mean of the first two seconds of the
%    gyroscope signals.
mu1 = mean(gyro_thigh_y(1:2*fs));
mu2 = mean(gyro_shank_y(1:2*fs));

% 10) Initialise the state vector.
x = [0, -(l1+l2), -85, 0, 0, -5, 0, 0, mu1, mu2]';

% 11) Initialise the error covariance matrix.
P = diag(ones(1, 10) * 1);

% 12) Define the measurement matrix.
H = [0 0 0 1 0 0 0 0 1 0; ...
     0 0 0 1 0 0 1 0 1 1; ...
     0 0 1 0 0 0 0 0 0 0; ...
     0 0 1 0 0 1 0 0 0 0];

% 13) Define process noise covariance matrix.
sigma_d = 1000;
sigma_t1 = 0.09;
sigma_t2 = 0.09;
sigma_b = 0.00000000001;
Q = [...
sigma_d 0 0           0             0      0 0 0 0 0; ...
0 sigma_d 0           0             0      0 0 0 0 0; ...
0 0 sigma_t1^9/9 sigma_t1^4/4 sigma_t1^5/5 0 0 0 0 0; ...
0 0 sigma_t1^4/4 sigma_t1^3/3 sigma_t1^2/2 0 0 0 0 0; ...
0 0 sigma_t1^5/5 sigma_t1^2/2   sigma_t1   0 0 0 0 0; ...
0 0 0 0 0 sigma_t2^9/9 sigma_t2^4/4 sigma_t2^5/5 0 0; ...
0 0 0 0 0 sigma_t2^4/4 sigma_t2^3/3 sigma_t2^2/2 0 0; ...
0 0 0 0 0 sigma_t2^5/5 sigma_t2^2/2   sigma_t2   0 0; ...
0 0 0 0 0       0             0          0 sigma_b 0; ...
0 0 0 0 0       0             0          0 0 sigma_b];

% 14) Compute sample variance of the first two
%     seconds of the gyroscope signals.
sigma_1 = var(gyro_thigh_y(1:2*fs));
sigma_2 = var(gyro_shank_y(1:2*fs));

% 15) Define measurement noise covariance matrix.
sigma_f_1 = 4;
sigma_f_2 = 4;
sigma_s_1 = 2;
sigma_s_2 = 2;
R = [sigma_1    0       0        0; ...
        0    sigma_2    0        0; ...
        0       0    sigma_s_1   0; ...
        0       0       0    sigma_s_2];
    
% 16) Define matrix function f.
function f_k = f
    
	f_k = [- l1 * c_w * x(4) * sind(x(3)) ...
           - l2 * c_w * (x(4) + x(7)) ...
           * sind(x(3) + x(6)); ...
           - l1 * c_w * x(4) * cosd(x(3)) ...
           - l2 * c_w * (x(4) + x(7)) ...
           * cosd(x(3) + x(6)); ...
           x(4);
           x(5);
           0;
           x(7);
           x(8);
           0;
           0;
           0];

end

% 17) Define Jacobian of F.
function F_k = F

    A = - l1 * c_w * x(4) * cosd(x(3)) ...
        - l2 * c_w * (x(4) + x(7)) * cosd(x(3) + x(6));
    B = + l1 * c_w * x(4) * sind(x(3)) ...
        + l2 * c_w * (x(4) + x(7)) * sind(x(3) + x(6));
    C = - l1 * sind(x(3)) - l2 * sind(x(3) + x(6));
    D = - l1 * cosd(x(3)) - l2 * cosd(x(3) + x(6));
    E = - l2 * c_w * (x(4) + x(7)) * cosd(x(3) + x(6));
    F = + l2 * c_w * (x(4) + x(7)) * sind(x(3) + x(6));
    G = - l2 * sind(x(3) + x(6));
    h = - l2 * cosd(x(3) + x(6));

    F_k = [0, 0, A, C, 0, E, G, 0, 0, 0; ...
           0, 0, B, D, 0, F, h, 0, 0, 0; ...
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0; ...
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0; ...
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0; ...
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0; ...
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...  
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

end

% 18) Filter loop.
for i=1:1:len
    
    % TIME UPDATE %
    
    % Compute fundamental matrix.
    Phi = eye(10) + F * Ts;
    
    % Compute a priori state estimate.
    x = x + f * Ts;
    
    % Compute a priori error covariance matrix.
    P = Phi * P * Phi' + Q;
    
    % CORRECT SENSOR READINGS %
    
    % Compute acceleration in the world frame that
    % occurs due to motion, based on a priori state
    % estimate.
    ax = - l1 * ((c_w * x(4))^2 * cosd(x(3)) ...
         + c_w * x(5) * sind(x(3))) ...
         - l2 * ((c_w * x(4) + c_w * x(7))^2 ...
         * cosd(x(3) + x(6)) ...
         + c_w * (x(5) + x(8)) * sind(x(3) + x(6)));
    az = - l1 * (c_w * x(5) * cosd(x(3)) ...
         - (c_w * x(4))^2 * sind(x(3))) ...
         - l2 * (c_w * (x(5) + x(8)) ...
         * cosd(x(3) + x(6)) ...
         - (c_w * x(4) + c_w * x(7))^2 ...
         * sind(x(3) + x(6)));
    
    % Normalise acceleration to gravity.
    ax_n = ax / gravity;
    az_n = az / gravity;
    
    % Compute transformation matrix
    Tz = [cosd(x(3) + x(6) + 90), 0, ...
         -sind(x(3) + x(6) + 90); 0, 1, 0; ...
          sind(x(3) + x(6) + 90), 0, ...
          cosd(x(3) + x(6) + 90)];
     
    % Rotate acceleration to body frame.
    a_mb = Tz * [ax_n; 0; az_n];

    % Compute gravity estimate by subtracting 
    % motion-based acceleration from sensor readings.
    g_c = [acc_shank_x(i); 0; acc_shank_z(i)] - a_mb;
     
    % Constitute the measurement vector from the
    % gyroscope signals, theta_1, and the corrected 
    % angle estimate theta_1 + theta_2.
    z = [gyro_thigh_y(i); gyro_shank_y(i); 0; 0];
    z(3) = atan2d(acc_thigh_z(i), acc_thigh_x(i)) - 180;
    z(4) = atan2d(g_c(3), g_c(1)) - 180;
    
    % Output map.
    theta12_c(i) = z(4);
    a_m(:, i) = a_mb;
    
    % Set sigma_3 in measurement noise covariance
    % matrix according to motion intensity.
    if marker(i) == 1
       R(3, 3) = sigma_f_1;
       R(4, 4) = sigma_f_2;
    end
    
    if marker(i) == 0
       R(3, 3) = sigma_s_1;
       R(4, 4) = sigma_s_2;
    end
    
    % MEASUREMENT UPDATE %
    
    % Compute Kalman gain.
    K = P * H' / (H * P * H' + R);
    
    % Compute a posteriori estimate.
    x = x + K * (z - H * x);
    
    % Update error covariance matrix.
    P = (eye(10) - K * H) * P;
    
    % Map internal states to output vector.
    theta1(i) = x(3);
    theta2(i) = x(6);
    
    X(:, i) = x;
    
end

end
