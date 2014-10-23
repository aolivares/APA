function F = eofKalman(p)

% FUNCTION EOFKALMAN Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the Kalman Filter. 

% *************************************************************************
% - Author: Alberto Olivares-Vicente.
% - Entity: University of Granada, Spain. 
% - Last revision: 08/21/2013.
% *************************************************************************
gw = gwLibrary;
% -------------------------------------------------------------------------
% 1) Set variables.
% -------------------------------------------------------------------------
global gyro;
global acc;
global frec;
global true_angle;

alpha = p(1);
beta = p(2);

angle_KF = gw.fusion_KF(gyro,acc,frec,var(acc),var(acc),var(gyro),alpha,...
    beta,acc(1));

F = sqrt(mean((true_angle-angle_KF').^2));
%--------------------------------------------------------------------------
% End of file \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%--------------------------------------------------------------------------