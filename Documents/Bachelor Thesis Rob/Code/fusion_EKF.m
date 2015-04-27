function [theta1, theta2] = fusion_EKF(gyro_thigh, ...
          gyro_shank, acc_thigh, acc_shank, fs, l1, l2)

% FUNCTION FUSION_EKF Applies an extended Kalman filter
% in order to fuse the accelerometer and gyroscope data
% and thus obtain an accurate orientation estimate of 
% the thighs and shanks.
%
% Input parameters:
% |_ 'gyro_thigh':  Vector containing the gyroscope 
%                   signal of the thigh.
% |_ 'gyro_shank':  Vector containing the gyroscope
%                   signal of the shank.
% |_ 'acc_thigh':   Vector containing the acceleration 
%                   signal of the thigh.
% |_ 'acc_shank':   Vector containing the acceleration
%                   signal of the shank.
% |_ 'frec':        Sampling frecuency. Must be real
%                   positive.
% |_ 'l1':          Length of thigh. Must be real
%                   positive.
% |_ 'l2':          Length of shank. Must be real
%                   positive.
%
% Output parameters:
% |_ 'theta1':      Thigh angle with respect to the
%                   x-axis of the world frame.
% |_ 'theta2':      Shank angle with respect to the
%                   thigh.
%
% IMPORTANT NOTE:   gyro_thigh, gyro_shank, acc_thigh,
%                   and acc_shank must have the same 
%                   length. Otherwise, an error will
%                   be returned.
% -----------------------------------------------------
% Authors: Robin Weiss
% Entity: University of Applied Sciences Munster,
%         Munster, Germany.
% Last modification: 26/04/2015.
% -----------------------------------------------------

% 1) Check input arguments.
if ~(gyro_thigh == gyro_shank == acc_thigh == acc_shank)    
    error(['Input arguments ''gyro_thigh'', ', ... 
           '''gyro_shank'', ''acc_thigh'', and ' ...
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

% 1) Import GaitWatch functions library. All existing
%    functions have to be called using 'gw.functionName'.
gw = gwLibrary;

% 2) Define the sampling period.
Ts = 1 / fs;

% Compute length of the signals.
len = length(gyro_thigh);

% Define rotation matrix.
syms Tz(theta)
Tz(theta) = [cos(theta),  0, sin(theta); ...
                  0,      1,      0;     ...
             -sin(theta), 0; cos(theta)];
         
% Compute intensity level.
lwin_fsd = 20;    
threshold_fsd = 3;    
shift_fsd = 19;    
lambda = 30;
[V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd, ...
                       shift_fsd, 512, threshold_fsd);
[marker,T_fsd_expanded] = gw.compEstMark(V_fsd, ...
                T_fsd,input_signal,lwin_fsd,shift_fsd);

% INITIALISATION OF PARAMETERS %

% 3) Compute mean of gyroscope signals of the first two
%    seconds.
mu1 = mean(gyro_thigh(1:2*fs));
mu2 = mean(gyro_thigh(1:2 * fs));

% 4) Initialise the state vector.
x = [0, -(l1+l2), -pi/2, 0, 0, 0, 0, 0, 0, mu1, mu2]';

% Map gyroscope signals to measurement vector. 
z = [gyro_thigh; gyro_shank; zeros(1, len)];

% 5) Initialise the covariance matrix.
P = diag(ones(1, 10));

% 7) Define the measurement matrix.
H = [0 0 0 1 0 0 0 0 1 0; ...
     0 0 0 1 0 0 1 0 0 1; ...
     0 0 1 0 0 1 0 0 0 0];

% 6) Define process noise covariance matrix.
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

% 8) Define measurement noise covariance matrix.
sigma_f = 0;
sigma_s = 0;
R = [1 0 0; ...
     0 1 0; ...
     0 0 sigma_s];

% 8) Define matrix function f.
syms f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = ...
[-l1 * x4 * sin(x3) - l2 * (x4 + x7) * sin(x3 + x6); ...
 -l1 * x4 * cos(x3) - l2 * (x4 + x7) * cos(x3 + x6); ...
 x4; x6; 0; x7; x8; 0; 0; 0];
         
% 7) Define matrix function F.
syms F(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
F(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = ...
[0, 0, 0, -l2 * sin(x3) - l2 * sin(x3+x6), 0, 0, ...
 - l2 * sin(x3+x6); ...
 0, 0, 0, l1 * cos(x3) + l2 * cos(x3+x6), 0, 0, ...
 + l2 * cos(x3+x6); ...
 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; ...
 0, 0, 0, 0, 1, 0, 0, 0, 0, 0; ...
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; ...
 0, 0, 0, 0, 0, 0, 0, 1, 0, 0; ...
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% Initialise output vectors
theta1 = zeros(1, len);
theta2 = zeros(1, len);

% 7) Filter loop.
for i=1:1:len-1
    
    % Set sigma3 in measurement noise covariance matrix
    % according to motion intensity.
    if marker(i)==1
       R(3, 3) = sigma_f;
    end
    if marker(i)==0
       R(3, 3) = sigma_s;
    end
    
    % TIME UPDATE %
    
    % Compute fundamental matrix.
    Phi = eye(10) * F(x) * Ts;
    
    % Compute a priori state estimate.
    x = x + f(x) * Ts;
    
    % Compute a priori error covariance matrix.
    P = Phi * P * Phi' + Q;
    
    % CORRECT SENSOR READINGS %
    
    % Compute acceleration due to motion.
    ax = -l1 (x(4)^2 * cos(x(3)) + x(5) sin(x(3)) ...
         - l2 ((x(4) + x(7))^2 * cos(x(3) + x(6)) ...
         + (x(5) + x(8)) sin(x(3) + x(6)));
    az = -l1 (x(5) * cos(x(3)) - x(4)^2 sin(x(3)) ...
         - l2 (x(5) + x(8)) * cos(x(3) + x(6)) ...
         + ((x(4) + x(7))^2 sin(x(3) + x(6)));
     
    % Rotate acceleration to body frame.
    a_rad = Tz(theta1 + theta2 - 2 * pi) * ax;
    a_tan = Tz(theta1 + theta2 - 2 * pi) * az;

    % Compute gravity estimate.
    g = acc_shank(i) - [a_rad; a_tan];
    
    % Compute corrected angle estimate and update 
    % measurement vector.
    z(3, i) = atan2(g(1), g(2));
    
    % MEASUREMENT UPDATE %
    
    % Compute Kalman gain.
    K = P * H' / (H * P * H' + R);
    
    % Compute a posteriori estimate.
    x = x + K * (z - H * x);
    
    % Update error covariance matrix.
    P = (eye(10) - K * H) / P;
    
    % Map internal states to output vector.
    theta1(i) = x(3);
    theta2(i) = x(6);
    
end

end