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
global true_angle;
global Ini;

gyro = gyro_KF;
acc = accM_angle_KF;
true_angle = reference;
frec = f;
Ini = accM_angle_KF(1);

% -------------------------------------------------------------------------
% 2) Call the minimization routine.
% -------------------------------------------------------------------------

disp('Optimizing parameters of Kalman Filter (it may take a while) ...');

xinit = p0_KF;
tol= 10^-6;
max_feval = 5000;
[xmin,fmin,ct] = ANMS(@eofKalman, xinit, tol, max_feval);
fprintf('\nOPTIMIZATION RESULTS\n--------------------\n')
fprintf('RMSE:%0.5f\nOPTIMAL PARAMETERS:\n - Alpha: %0.9f\n - Beta: %0.9f\n',...
    fmin,xmin(1),xmin(2))
%--------------------------------------------------------------------------
% End of file \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%--------------------------------------------------------------------------