function [theta1, theta2] = fusion_EKF(gyro_thigh, ...
          gyro_shank, acc_thigh, acc_shank, fs, a1, a2)

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
% |_ 'a1':          Length of thigh.
% |_ 'a2':          Length of shank.
%
% Output parameters:
% |_ 'theta1':      Thigh angle.
% |_ 'theta2':      Shank angle with respect to the
%                   thigh.
%
% -----------------------------------------------------
% Authors: Robin Weiss
% Entity: University of Applied Sciences Munster
% Last modification: 26/04/2015.
% -----------------------------------------------------
len = length(gyro_thigh);

% 1) Define the sampling period.
Ts = 1 / fs;

% 2) Compute mean of gyroscope signals of the first two
%    seconds.
mu1 = mean(gyro_thigh(1:2*fs));
mu1 = mean(gyro_thigh(1:2 * fs));

% 3) Initialise the state vector.
x = [0, -(a1+a2), -pi/2, 0, 0, 0, 0, 0, 0, mu1, mu2]';

% 2) Initialise the covariance matrix.
P = diag(ones(1, 10));

% 3) Set process noise covariance matrix.
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

% 4) Define the measurement matrix.
H = [0 0 0 1 0 0 0 0 1 0; ...
     0 0 0 1 0 0 1 0 0 1; ...
     0 0 1 0 0 1 0 0 0 0];

% 6) Define function f.
syms f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = ...
[-a1 * x4 * sin(x3) - a2 * (x4 + x7) * sin(x3 + x6); ...
 -a1 * x4 * cos(x3) - a2 * (x4 + x7) * cos(x3 + x6); ...
 x4; x6; 0; x7; x8; 0; 0; 0];
         
% 7) Define matrix F.
F = jacobian(f, ...
            [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]');

% 7) Filter loop.
for i=1:1:len-1
    
    % Compute fundamental matrix.
    Phi = eye(10) * F(x) * Ts;
    
    % Compute a priori state estimate.
    x = x + f(x) * Ts;
    
end

end