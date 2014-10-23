% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||| KALMAN FILTER PARAMETER OPTIMIZATION |||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% The following file contains the routine which uses an orientation
% reference system (in this case a Qualisys optical motion tracking system)
% to optimize the parameters of the orientation estimation algorithms
% applied on the signals gathered by the GaitWatch system. 
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai Bötzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universität München.
% Last modification: 13/12/2013.
% -------------------------------------------------------------------------

clear all;close all;clc;

gw = gwLibrary;

% 1) Load Qualisys and GaitWatch data and synchronize them.
% -------------------------------------------------------------------------
% 1.1) We first call the main routine to load and extract Qualisys data.
readQualysis;

% 1.2) We then call the main routine to load, calibrate and analyze 
%      GaitWatch data.
main;

% 1.3) Select the GW signal that is going to be used as a reference to
%      synchronize GW data with QS data.
integ_gyro_variables = who('-var', 'pitch_gyro*');
[Selection2] = listdlg('ListString',integ_gyro_variables,'Name',...
        'Select the GW synchronization signal','ListSize',[300 120],...
        'SelectionMode','single');
gw_signal = eval(integ_gyro_variables{Selection2});

% 1.4) Select the QS signal that is going to be used as a reference to
%      synchronize GW data with QS data.
qs_angle_variables = who('-var', '*_qs');
[Selection3] = listdlg('ListString',qs_angle_variables,'Name',...
        'Select the QS synchronization signal','ListSize',[300 120],...
        'SelectionMode','single');
qs_signal = eval(qs_angle_variables{Selection3});

% 1.5) Plot the two figures so the user can select the synchronization
%      indexes.
[index1, index2] = gw.sync_GW_QS(gw_signal, qs_signal);

% There are two different kinds of body segments which parameters can be
% optimized: the leg sements and the trunk.

% 2) If the user has selected leg segments.
% -------------------------------------------------------------------------
if sum(strcmpi(S(Selection),'trunk')) == 0
    
    % 2.1) Synchronize signals.
    gw_signal = gw_signal(index1 - index2: index1 - index2 + ...
        length(qs_signal)-1);
    ax = ax(index1 - index2: index1 - index2 + length(qs_signal)-1);
    az = az(index1 - index2: index1 - index2 + length(qs_signal)-1);
    gy = gy(index1 - index2: index1 - index2 + length(qs_signal)-1);
    
    % 2.2) If the Qualisys roll signal exists, then we need to compensate
    %      the acceleration in the Z axis prior to the computation of the
    %      pitch using the GaitWatch data.
    if exist('roll_qs','var')
        
        figure
        plot(time_qs,az)
        hold on
        plot(time_qs,cosd(roll_qs).*az','--r')
        legend('Uncompensated','Compensated')

        % Recompute pitch with roll-compensated acceleration.
        pitch_acc = atan2d(cosd(roll_qs).*az',ax');
        pitch_acc = gw.correct_quad_shifts(pitch_acc,'deg') - 90;
    
    % 2.3) If the Qualisys signal does not exist, then the pitch angle is
    %      directly computed without any roll compensation.
    elseif ~exist('roll_qs','var')
        pitch_acc = atan2d(az,ax);
        pitch_acc = gw.correct_quad_shifts(pitch_acc,'deg') - 90;
    end
    % 2.4) Plot synchronized signals (Qalisys pitch and GaitWatch pitch).
    figure
    plot(time_qs,pitch_acc)
    hold on
    plot(time_qs,qs_signal+(pitch_acc(1)-qs_signal(1)),'r')
    title('SYNCHRONIZED SIGNALS')
    legend('Gaitwatch','Qualisys')
    
    % 3) Optimize the Kalman filter parameters.
    % ---------------------------------------------------------------------
    
    reference = qs_signal+(pitch_acc(1)-qs_signal(1));
    accM_angle_KF = pitch_acc;
    accM_angle_GKF = pitch_acc;
    
    lwin_fsd = 20;    
    threshold_fsd = 3;    
    shift_fsd = 19;    
    lambda = 30;
    input_signal = sqrt(ax.^2+az.^2);
    [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,threshold_fsd);
    [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,input_signal,...
        lwin_fsd,shift_fsd);
            
    gyro_KF = gy;
    gyro_GKF = gy;

    % 3.1)Set initial value of parameters.
    p0_KF = [10000 10];
    p0_GKF = [100 1000 0.1 0.1];
    
    % 3.2) Call the KF optimization routine.
    optimizeKF;
    
    % 3.3) Extract optimal parameters.
    opt_alpha_KF = xmin(1);
    opt_beta_KF = xmin(2);

    % 3.4) Call the GKF optimization routine.
    optimizeGKF;
    
    % 3.5) Extract optimal parameters.
    opt_alpha1_GKF = xmin(1);
    opt_alpha2_GKF = xmin(2);
    opt_beta1_GKF = xmin(3);
    opt_beta2_GKF = xmin(4);
        
    % 3.4) Apply Kalman filter using optimal parameters.
    pitch_KF_right_shank = gw.fusion_KF(gyro_KF,accM_angle_KF,f,...
        var(accM_angle_KF),var(accM_angle_KF),...
        var(gyro_KF),opt_alpha_KF,opt_beta_KF,accM_angle_KF(1));

    % 3.5) Compute RMSE.
    rmse_KF = sqrt(mean((reference(1:end-1) - pitch_KF_right_shank(1:end-1)').^2));

    % 3.6) Plot results.
    figure
    hold on
    plot(time_qs,reference,'black')
    plot(time_qs,accM_angle_KF)
    plot(time_qs,pitch_KF_right_shank,'red')
    legend('Qualisys reference','Accelerometer',...
        sprintf('GaitWatch - Kalman Filter, Alpha=%0.4f, Beta=%0.4e, RMSE=%0.4f',...
        opt_alpha_KF,opt_beta_KF,rmse_KF))
    set(gca,'fontsize',14)
    xlabel('Time(s)','fontsize',18)
    ylabel('Angle (deg)','fontsize',18)
end           
        