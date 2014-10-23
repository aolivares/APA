% The goal of this routine is to use the Gauss-Newton Algorithm to minimize
% the error between the actual orientation angle and the one estimated with
% the Kalman Filter. 
% 
% *************************************************************************
% - Author: Alberto Olivares-Vicente.
% - Entity: University of Granada, Spain. 
% - Last revision: 08/21/2013.
% *************************************************************************

% -------------------------------------------------------------------------
% 1) Set variables.
% -------------------------------------------------------------------------
global gyro;
global acc;
global frec;
global pObsVar;
global accVar;
global gyroVar;
global true_angle;
global marker;

gyro = gyro_GKF;
acc = accM_angle_GKF;
pObsVar = var(accM_angle_GKF);
accVar = var(accM_angle_GKF);
gyroVar = var(gyro_GKF);
true_angle = reference;
frec = f;
marker = marker_fsd;

% -------------------------------------------------------------------------
% 2) Call the minimization routine.
% -------------------------------------------------------------------------
%clc;
disp('Optimizing parameters of Gated Kalman Filter (it may take a while) ...');

xinit = p0_GKF;
tol = 10^-6;
max_feval= 5000;
[xmin,fmin,ct] = ANMS(@eofGKF,xinit,tol,max_feval);
% [p,a]=lsqnonlin(@eofKalman,p0_KF3,-inf,inf,optimset('maxf',500000,...
%     'TolFun',1e-8,'TolX',1e-8,'MaxIter',1e4));
fprintf('\nOPTIMIZATION RESULTS\n--------------------\n')
fprintf('RMSE:%0.5f\nOPTIMAL PARAMETERS:\n - Alpha1: %0.9f\n - Alpha2: %0.9f\n - Beta1: %0.9f\n - Beta2: %0.9f\n',...
    fmin,xmin(1),xmin(2),xmin(3),xmin(4))
%--------------------------------------------------------------------------
% End of file \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%--------------------------------------------------------------------------