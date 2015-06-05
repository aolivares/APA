% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------- Gait Watch VS Qualisys System -----------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Veronica  Torres (2): 
%                   |_ vts24@correo.ugr.es 
%                 - Dr. Eng. Alberto Olivares (3): 
%                   |_ aolivares@ugr.es
%                 - Robin Weiss (4): 
%                   |_ mail@robinweiss.de
%
% * Affiliation: (1) Laboratory of Motion Analysis and Deep Brain 
%                    Stimulation, Department of Neurology, Klinikum 
%                    Grosshadern, Munich, Germany.
%                (2) Master in Telecommunication Engineering, University of 
%                    Granada, Granada, Spain, (student).
%                (3) Signal Processing and Biomedical Applications group,
%                    Department of Signal Theory, Telematics and
%                    Communications, University of Granada, Granada, Spain.
%                (4) Bachelor in Electrical Engineering, University of 
%                    Applied Sciences of Munster, Munster, Germany, 
%                    (student).
%
% * Last modification: 27/05/2015

% INFORMATION: This file contains the routine to calculate the orientation
% angle in the GW and QS and the comparation between both.
% This script is  a modified version of the 'main.m' and 'main_analysis.m'
% file implemented by Alberto Olivares.
% 
% * 1) Extract data from both systems.
% 
% * 2) Calculate pitch in GW System.
% 
% * 3) Calculate pitch in QS System.
% 
% * 4) Comparation.
%
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

comparation = 'no';
showPlots = 'yes';

% Suppress warnings if no peak is detected during the calibration.
  warning('off', 'signal:findpeaks:largeMinPeakHeight')
  
% -------------------------------------------------------------------------
% 1) Extract data.
% -------------------------------------------------------------------------

% Import GaitWatch functions library.
gw = gwLibrary;

% -------------------------------------------------------------------------
% 1.2) Define structure of data array.
% -------------------------------------------------------------------------
% This array determines the behavior of the all the code which is
% subsequently executed. It contains the data format of the file stored in
% the GaitWatch memory card. The first element indicates the channel 
% number, the second element indicates the measured magnitude ('a' for
% acceleration, 'g' for angular rate and 'h' for magnetic field), the third
% element identifies the position of the body segment (either right or 
% left, or center, in case the sensors are placed close to the intersection
% of the Sagittal and Coronal Planes), the fifth element indicates the body
% segment to which the sensor is attached, and finally, the sixth element
% indicates the kind of calibration process which is applied to the raw
% data (1: uniaxial, 3: triaxial).

data_struct={1, 'g','Y',  'right',  'shank', 1;
             2, 'a','Z',  'right',  'shank', 1;
             3, 'a','X',  'right',  'shank', 1;
             4, 'g','Y',  'right',  'thigh', 1;
             5, 'a','Z',  'right',  'thigh', 1;
             6, 'a','X',  'right',  'thigh', 1;
             7, 'g','Y',  'left',   'shank', 1;
             8, 'a','Z',  'left',   'shank', 1;
             9, 'a','X',  'left',   'shank', 1;
             10,'g','Y',  'left',   'thigh', 1;
             11,'a','Z',  'left',   'thigh', 1;
             12,'a','X',  'left',   'thigh', 1;
             13,'g','Y',  'left',   'arm'  , 1;
             14,'g','X',  'left',   'arm'  , 1; 
             15,'g','Y',  'right',  'arm'  , 1;
             16,'g','X',  'right',  'arm'  , 1; 
             17,'g','Z',  'center', 'trunk', 1;
             18,'g','Y',  'center', 'trunk', 1;
             19,'g','X',  'center', 'trunk', 1;
             20,'a','Z',  'center', 'trunk', 3;
             21,'a','Y',  'center', 'trunk', 3;
             22,'a','X',  'center', 'trunk', 3;
             23,'h','X',  'center', 'trunk', 3;
             24,'h','Y',  'center', 'trunk', 3;
             25,'h','Z',  'center', 'trunk', 3;
             };
         
size_data_struct=size(data_struct);    

% -------------------------------------------------------------------------
% 2.2) Load data file.
% -------------------------------------------------------------------------

% Select and load data from the hard drive.
[filename_GW, filepath] = uigetfile('*.mat', ...
    'Select the GW data file (.mat)', '../../../Treadmill experiments/GW data/last recording treadmill');
load(fullfile(filepath ,char(filename_GW)));

% Build vector containing time samples.
[f, date, start_time, end_time, file_id] = gw.getFHinfo(FileHeader);

% Build time signal. 
len_data = length(data);
time = (0:len_data-1) / f;

% Reshape data (split channel 23 into 3 channels and append them to the
% data matrix)
[mag_x, mag_y, mag_z] = gw.getMagData(double(data(:,23)));

% Compute the sampling frequency of the magnetometer.
f_mag = f/(length(data)/length(mag_x));

% We now build a time signal using the newly computed frequency. This
% signal will only be used in figures (to plot magnetic field vs. time). 
time_mag=zeros(1,length(mag_x));
for i=1:length(mag_x)-1
 time_mag(i+1)=time_mag(i)+1/f_mag;
end

% Correct magnetometer signals.
% mag_x = gw.correct_mag_data(mag_x);
% mag_y = gw.correct_mag_data(mag_y);
% mag_z = gw.correct_mag_data(mag_z);
% close all;

% -------------------------------------------------------------------------
% 2.3) Correct magnetometer signals automatically.
% -------------------------------------------------------------------------
% Magnetometer X:

% Detect the noise peaks using a threshold far above the mean. We detect
% the positives peaks as well as negatives.
mean_mag_x = mean(mag_x);
threshold_x = 2000 + mean_mag_x;

[~ , locs_high] = findpeaks(mag_x, 'minpeakheight', threshold_x);
[~ , locs_short] = findpeaks(- mag_x, 'minpeakheight', threshold_x);

indexes  =[locs_high locs_short];

% Substitute the erroenous values by the subsequent value (unless the
% erroneous value is in the last position of the singal, in which case 
% it is substituted by the preceding value). 

for i = 1:length(indexes)
    if indexes(i) == length(mag_x)
        mag_x(indexes(i)) = mag_x(indexes(i)-1);
    else
        mag_x(indexes(i)) = mag_x(indexes(i)+1);
    end
end

% Magnetometer Y:

% Detect the noise peaks using a threshold far above the mean.
mean_mag_y = mean(mag_y);
threshold_y = 1000 + mean_mag_y;

[~ , locs_high] = findpeaks(mag_y, 'minpeakheight', threshold_y);
[~ , locs_short] = findpeaks(- mag_y, 'minpeakheight', threshold_y);

indexes  =[locs_high locs_short];

% Substitute the erroenous values by the subsequent value (unless the
% erroneous value is in the last position of the singal, in which case 
% it is substituted by the preceding value). 

for i = 1:length(indexes)
    if indexes(i) == length(mag_y)
        mag_y(indexes(i)) = mag_y(indexes(i)-1);
    else
        mag_y(indexes(i)) = mag_y(indexes(i)+1);
    end
end

% Magnetometer Z:

% Detect the noise peaks using a threshold far above the mean.
mean_mag_z = mean(mag_z);
threshold_z = 1000 + mean_mag_z;

[~ , locs_high] = findpeaks(mag_z, 'minpeakheight', threshold_z);
[~ , locs_short] = findpeaks(- mag_z, 'minpeakheight', threshold_z);

indexes  =[locs_high locs_short];

% Substitute the erroenous values by the subsequent value (unless the
% erroneous value is in the last position of the singal, in which case 
% it is substituted by the preceding value). 

for i = 1:length(indexes)
    if indexes(i) == length(mag_z)
        mag_z(indexes(i)) =mag_z(indexes(i)-1);
    else
        mag_z(indexes(i)) = mag_z(indexes(i)+1);
    end
end

% Interpolate the magnetometer signals
mag_x_interp = interp1(time_mag,mag_x,time,'spline');
mag_y_interp = interp1(time_mag,mag_y,time,'spline');
mag_z_interp = interp1(time_mag,mag_z,time,'spline');

% if strcmpi(showMagData,'yes')
%     figure
%     subplot(3,1,1)
%     plot(time_mag,mag_x)
%     hold on
%     plot(time,mag_x_interp,'--r')
%     legend('Original','Interpolated')
%     xlabel('Time (s)');
%     ylabel('Magnetic field (raw)')
%     subplot(3,1,2)
%     plot(time_mag,mag_y)
%     hold on
%     plot(time,mag_y_interp,'--r')
%     legend('Original','Interpolated')
%     xlabel('Time (s)');
%     ylabel('Magnetic field (raw)')
%     subplot(3,1,3)
%     plot(time_mag,mag_z)
%     hold on
%     plot(time,mag_z_interp,'--r')
%     legend('Original','Interpolated')
%     xlabel('Time (s)');
%     ylabel('Magnetic field (raw)')
% end

% Remove channel 23 from the data matrix and add three new channels 
% containing the interpolated magnetometer signals to the data matrix.
data = [double(data(:,1:22)) mag_x_interp' mag_y_interp' mag_z_interp'];

% -------------------------------------------------------------------------
% 2.4) Identify, extract and plot data.
% -------------------------------------------------------------------------
% We now have to extract the data from the 'data' matrix. This is done 
% dynamically, so the code is flexible and only depends on 'data_struct'. 
% A vector is created for each one of the magnitudes which were measured by
% the unit placed on the segment selected by the user from the list.

% Show segment selection menu (multiple selection by holding the CTRL key
% is also possible). 
select_ok_flag = 0;
while select_ok_flag == 0
    
    % Define the list of segments which are shown to the user. 
    S = cell(1,12);
    for i = 1:12 % Only Select shanks and thighs.
        S{i} = [data_struct{i,4},' ',data_struct{i,5}];
    end
    S = unique(S,'stable');
    
    % Show the selection dialog to the user.
    [Selection,ok] = listdlg('ListString',S,'Name',...
        'Select the unit you wish to calibrate','ListSize',[160 120],...
        'SelectionMode','multiple');
    if ~isempty(Selection) 
        select_ok_flag = 1;
    else
        msg = msgbox('Please select at least one segment');
        uiwait(msg);
    end
end

% We now extract all the data channels and dinamically generate the
% variable name according to the information in 'data_struct' and the
% selection made by the user.

size_data_struct=size(data_struct);

for i = 1:length(Selection)
    
    % We split the selection in a two-elements cell.
    position_segment = strsplit(S{Selection(i)});
    position = position_segment{1};
    segment = position_segment{2};
    
    for j=1:size_data_struct(1)
        
        % We now want to extract any data which were measured for the
        % selected segment(s) and position(s).
        if strcmpi(data_struct{j,4},position) && ...
                strcmpi(data_struct{j,5},segment)
                
                % We build the name of the variable following this format:
                % 'magnitude_axis_position_segment_calibrationType'
                var_name=strcat(data_struct{j,2},'_',data_struct{j,3},...
                    '_',data_struct{j,4},'_',data_struct{j,5},'_',...
                num2str(data_struct{j,6}));
                
                % We build the code that we wish to execute.
                code_line=strcat(var_name,'=data(:,%d);');
                
                % We evaluate the code line.
                eval(sprintf(code_line,j));
        end
    end
end

% -------------------------------------------------------------------------
% 2.5) Calibrate acceleration.
% -------------------------------------------------------------------------
% We now search all the variables containing raw acceleration. To do so, we
% build a list with all variables starting by 'a_'.
% IMPORTANT: ONLY VARIABLES CONTAINING RAW ACCELERATION DATA CAN HAVE A 
% NAME STARTING BY 'a_', OTHERWISE THE ROUTINE WILL FAIL. 

acc_variables = who('-var', 'a_*');
acc_variables_3 = {};
if ~isempty(acc_variables)
    
    for i = 1: length(acc_variables)
        var_name = acc_variables{i};

        % Uniaxial calibration.
        if strcmp(var_name(end),'1')
            split_var_name = strsplit(var_name,'_');

            var_values = eval(var_name);
            var_name_axis = split_var_name{2};
            var_name_pos = split_var_name{3};
            var_name_seg = split_var_name{4};
            cal_var_name = strcat(var_name,'_C');
            code_line = strcat(cal_var_name,['= gw.acc1DCal(var_values',...
                ',var_name_axis,var_name_pos,var_name_seg);']);
            eval(code_line);

        end

        % Store name of variables needing triaxial calibration.
        if strcmp(var_name(end),'3')
            acc_variables_3{i} = acc_variables{i};
        end

    end

    if ~isempty(acc_variables_3)

        % Remove empty elements from the 'acc_variables_3' list.
        acc_variables_3 = acc_variables_3(~cellfun('isempty',...
            acc_variables_3));

        % If there are more than 2 segments with a triaxial accelerometer, 
        % then the cell 'acc_variables_3' would look like: 
        % [a_X_pos_segment1_3, a_X_pos_segment2_3, a_Y_pos_segment1_3,
        %  a_Y_pos_segment2_3, a_Z_pos_segment1_3, a_Z_pos_segment2_3].
        % That is a (1 x 3N) vector where N is the number of segments.
        %
        % We actually want it to look like like:
        % [a_X_pos_segment1_3,  a_Y_pos_segment1_3, a_Z_pos_segment1_3;
        %  a_X_pos_segment2_3,  a_Y_pos_segment2_3, a_Z_pos_segment2_3], 
        % that is a (N x 3) matrix. 
        % This is done by reshaping the cell.

        % Check if there are more than 3 variables (which would mean that 
        % there are more than two segments with triaxial accelerometers). 
        % If so, reshape the list. 
        if length(acc_variables_3) > 3
            acc_variables_3 = reshape(acc_variables_3,...
                length(acc_variables_3)/3,3);
        end

        size_acc_variables_3 = size(acc_variables_3);
        
        % We now sweep every row of the list and build the new variables
        % which will contain the calibrated data.
        for i = 1:size_acc_variables_3(1)
            split_var_name = strsplit(acc_variables_3{i},'_');
            var_name_seg = split_var_name{4};
            
            % We append '_C' at the end of the existing variable names.
            cal_var_name_x = strcat(acc_variables_3{i,1},'_C');
            var_values_x = eval(acc_variables_3{i,1});
            cal_var_name_y = strcat(acc_variables_3{i,2},'_C');
            var_values_y = eval(acc_variables_3{i,2});
            cal_var_name_z = strcat(acc_variables_3{i,3},'_C');
            var_values_z = eval(acc_variables_3{i,3});
            
            % Build the code line which, when evaluated, calls the
            % calibration routine.
            code_line = strcat('[',cal_var_name_x,',',cal_var_name_y,',',...
                cal_var_name_z,[']=gw.acc3DCal(var_values_x,var_values_y,',...
                'var_values_z,var_name_seg);']);
            eval(code_line);
        end
    end
end

% -------------------------------------------------------------------------
% 2.6) Calibrate angular rate.
% -------------------------------------------------------------------------
% We now search all the variables containing raw angular rate. To do so, we
% build a list with all variables starting by 'g_'.
% IMPORTANT: ONLY VARIABLES CONTAINING RAW ANGULAR RATE DATA CAN HAVE A 
% NAME STARTING BY 'g_', OTHERWISE THE ROUTINE WILL FAIL. 

gyro_variables = who('-var', 'g_*');

if ~isempty(gyro_variables)
    for i = 1: length(gyro_variables)
        var_name = gyro_variables{i};

        % Uniaxial calibration.
        if strcmp(var_name(end),'1')
            
            % Split the name of the variable and get the axis, position and
            % segment.
            split_var_name = strsplit(var_name,'_');
            var_values = eval(var_name);
            var_name_axis = split_var_name{2};
            var_name_pos = split_var_name{3};
            var_name_seg = split_var_name{4};
            
            % Append '_C' to the variable name.
            cal_var_name = strcat(var_name,'_C');
            
            % Build the code line which, when evaluated, calls the
            % calibration routine.
            code_line = strcat(cal_var_name,['= gw.gyro1DCal(var_values,',...
                'var_name_axis,var_name_pos,var_name_seg);']);
            eval(code_line);
        end
    end
end


% -------------------------------------------------------------------------
% 2.7) Calibrate magnetic field.
% -------------------------------------------------------------------------
% We now search all the variables containing raw magnetic field. To do so, 
% we build a list with all variables starting by 'h_'.
% IMPORTANT: ONLY VARIABLES CONTAINING RAW MAGNETIC FIELD DATA CAN HAVE A 
% NAME STARTING BY 'h_', OTHERWISE THE ROUTINE WILL FAIL. 
mag_variables = who('-var', 'h_*')';

if ~isempty(mag_variables)
    % Again, if there are more than 3 variables (2 units) we reshape the 
    % list.
    if length(mag_variables) > 3
        mag_variables = reshape(mag_variables,length(mag_variables)/3,3);
    end

    size_mag_variables = size(mag_variables);
    for i = 1:size_mag_variables(1)
            split_var_name = strsplit(mag_variables{i},'_');
            var_name_seg = split_var_name{4};

            cal_var_name_x = strcat(mag_variables{i,1},'_C');
            var_values_x = eval(mag_variables{i,1});
            cal_var_name_y = strcat(mag_variables{i,2},'_C');
            var_values_y = eval(mag_variables{i,2});
            cal_var_name_z = strcat(mag_variables{i,3},'_C');
            var_values_z = eval(mag_variables{i,3});

            code_line = strcat('[',cal_var_name_x,',',cal_var_name_y,',',...
                cal_var_name_z,[']=gw.mag3DCal(var_values_x,var_values_y,',...
                'var_values_z,var_name_seg);']);
            eval(code_line);
    end
end


% -------------------------------------------------------------------------
% 2.8) Compute orientation.
% -------------------------------------------------------------------------
% Once we have extracted and calibrated all the data from the selected
% segments, we proceed to compute the orientation of the segment with
% respect to the Earth's surface. The procedure to computate the 
% orientation is different depending o the kind and amount of available
% sensors at the units placed on the segments. 
% As explained in the manual, we have three different kinds of units:
%
%  - Shank and thigh units: They contain a biaxial accelerometer (X and Z
%    axes) and a uniaxial gyroscope (Y axis). With these sensors we can
%    only compute the pitch angle. The computation of the algorithm is done
%    by fusing the angle estimation calculated using the accelerometers and
%    the angular rate measured by the gyroscope. The fusion is done using a
%    Kalman filter. 
%
%  - Arms: They only contain a biaxial gyroscope (X and Y axes). With this
%    sensor we can compute both pitch and roll angles. However, since we 
%    only have angular rate data, the only possible method to compute these
%    angles is to integrate the angular rate. Ideally this method would be
%    enough to obtain the roation angle, however, due to the presence of
%    noise in the output and the dynamic variation of the offset, the
%    computed angle drifts over time. A high pass filter can be applied to
%    partially remove the drift.
%
%  - Trunk: This unit contains a triaxial accelerometer, a triaxial
%    magnetometer and a triaxial gyroscope. Therefore, we can compute the
%    orientation with nine degrees of freedom (9DOF), that is, we can
%    estimate the pitch, roll and yaw angles. The estimation of these
%    angles is done by fusing the acceleration, magnetic field and angular
%    rate data using a quaternion extended kalman filter based on gradient 
%    descent.
    

for i = 1:length(Selection)
    
    position_segment = S{Selection(i)};
    switch position_segment
        case 'right shank'
            % Define the necessary signals.
            ax = a_X_right_shank_1_C;
            az = a_Z_right_shank_1_C;
            gy = g_Y_right_shank_1_C;
            
            % Compute pitch using acceleration.
            pitch_acc_right_shank = atan2d(az,ax);
            pitch_acc_right_shank = gw.correct_quad_shifts(...
                pitch_acc_right_shank,'deg') - 90;
            
            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_right_shank = gw.fusion_KF(gy, pitch_acc_right_shank,...
                f, var(pitch_acc_right_shank), var(pitch_acc_right_shank),...
                var(gy), alpha_KF, beta_KF, pitch_acc_right_shank(1));
            
            % Centre the signal
            pitch_KF_right_shank = pitch_KF_right_shank - pitch_KF_right_shank(1);

           % Calculte features to characterise the movement.
            [value_angle_neg,loc_angle_neg] = findpeaks(-pitch_KF_right_shank,'minpeakheight',7,'minpeakdistance',200);

            [value_angle_pos, loc_angle_pos] = findpeaks(pitch_KF_right_shank(...
                loc_angle_neg(1):length(pitch_KF_right_shank)), 'minpeakheight',7, 'minpeakdistance',200,...
                'Npeaks',length(loc_angle_neg));

            loc_angle_pos = loc_angle_pos + loc_angle_neg(1);
            stride_GW (1,1) = mean(diff(loc_angle_pos));
            stride_GW (2,1) = abs(mean(diff(loc_angle_pos)-stride_GW (1,1)));

%             swing_GW (1,1) = mean(loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1));
%             swing_GW (2,1) = abs(mean((loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1) )- swing_GW(1,1)));

            angle_GW(1,1) = mean(value_angle_pos);
            angle_GW(2,1) = mean(value_angle_neg);
            
            % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                input_signal,lwin_fsd,shift_fsd);
            
            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_right_shank = gw.fusionGKF(gy,pitch_acc_right_shank,...
                f,var(pitch_acc_right_shank),var(pitch_acc_right_shank),...
                var(gy),alpha1,alpha2,beta1,beta2,pitch_acc_right_shank(1),...
                marker_fsd);
            
            if strcmpi(comparation,'yes')
                % Integrate angular rate (just for comparation purposes).            
                ini_pos = pitch_acc_right_shank(1);
                pitch_gyro_right_shank = gw.integRate(1/f,gy,ini_pos);
                figure
                plot(time,pitch_acc_right_shank)
                hold on
                plot(time,pitch_gyro_right_shank,'r')
                plot(time,pitch_KF_right_shank ,'black')
                plot(time,pitch_GKF_right_shank,'cyan')
                title('Pitch angle of right shank - Comparison')
                xlabel('Time (s)')
                ylabel('Pitch (deg)')
                legend('Accelerometer-based','Angular rate integration',...
                    'Kalman filter','Gated Kalman Filter')              
%             else
%                 figure
%                 plot(time,pitch_KF_right_shank ,'black')
%                 title('Pitch angle of right shank')
%                 xlabel('Time (s)')
%                 ylabel('Pitch (deg)')
%                 legend('Kalman filter')
            end
            
        case 'left shank'
            % Define the necessary signals.
            ax = a_X_left_shank_1_C;
            az = a_Z_left_shank_1_C;
            gy = g_Y_left_shank_1_C;
            
            % Compute pitch using acceleration.
            pitch_acc_left_shank = atan2d(az,ax);
            pitch_acc_left_shank = gw.correct_quad_shifts(...
                pitch_acc_left_shank,'deg') - 90;
            
            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_left_shank = gw.fusion_KF(gy, pitch_acc_left_shank,...
                f, var(pitch_acc_left_shank), var(pitch_acc_left_shank),...
                var(gy), alpha_KF, beta_KF, pitch_acc_left_shank(1));
            
            % Centre the signal.
            pitch_KF_left_shank = pitch_KF_left_shank - pitch_KF_left_shank(1);
            
           % Calculte features to characterise the movement.
            [value_angle_neg,loc_angle_neg] = findpeaks(-pitch_KF_left_shank,'minpeakheight',7,'minpeakdistance',200);

            [value_angle_pos, loc_angle_pos] = findpeaks(pitch_KF_left_shank(...
                loc_angle_neg(1):length(pitch_KF_left_shank)), 'minpeakheight',7, 'minpeakdistance',200,...
                'Npeaks',length(loc_angle_neg));

            loc_angle_pos = loc_angle_pos + loc_angle_neg(1);
            stride_GW (1,3) = mean(diff(loc_angle_pos));
            stride_GW (2,3) =abs( mean(diff(loc_angle_pos)-stride_GW(1,3)));

%             swing_GW (1,2) = mean(loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1));
%             swing_GW (2,2) = abs( mean((loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1) )- swing_GW(1,2)));

            angle_GW(1,3) = mean(value_angle_pos);
            angle_GW(2,3) = mean(value_angle_neg);
            
             % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                input_signal,lwin_fsd,shift_fsd);
            
            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_left_shank = gw.fusionGKF(gy,pitch_acc_left_shank,...
                f,var(pitch_acc_left_shank),var(pitch_acc_left_shank),...
                var(gy),alpha1,alpha2,beta1,beta2,pitch_acc_left_shank(1),...
                marker_fsd);
            
            if strcmpi(comparation,'yes')
                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_left_shank(1);
                pitch_gyro_left_shank = gw.integRate(1/f,gy,ini_pos);
                figure
                plot(time,pitch_acc_left_shank)
                hold on
                plot(time,pitch_gyro_left_shank,'r')
                plot(time,pitch_KF_left_shank ,'black')
                plot(time,pitch_GKF_left_shank, 'cyan')
                title('Pitch angle of left shank - Comparison')
                xlabel('Time (s)')
                ylabel('Pitch (deg)')
                legend('Accelerometer-based','Angular rate integration',...
                    'Kalman filter','Gated Kalman filter')  
%             else
%                 figure
%                 plot(time,pitch_KF_left_shank ,'black')
%                 title('Pitch angle of left shank')
%                 xlabel('Time (s)')
%                 ylabel('Pitch (deg)')
%                 legend('Kalman filter')
            end
            
        case 'right thigh'
            % Define the necessary signals.
            ax = a_X_right_thigh_1_C;
            az = a_Z_right_thigh_1_C;
            gy = g_Y_right_thigh_1_C;
            
            % Compute pitch using acceleration.
            pitch_acc_right_thigh = atan2d(az,ax);
            pitch_acc_right_thigh = gw.correct_quad_shifts(...
                pitch_acc_right_thigh,'deg') - 90;
            
            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_right_thigh = gw.fusion_KF(gy, pitch_acc_right_thigh,...
                f, var(pitch_acc_right_thigh), var(pitch_acc_right_thigh),...
                var(gy), alpha_KF, beta_KF, pitch_acc_right_thigh(1));
            
            % Centre the signal.
            pitch_KF_right_thigh = pitch_KF_right_thigh - pitch_KF_right_thigh(1);
            
           % Calculte features to characterise the movement.
            [value_angle_neg,loc_angle_neg] = findpeaks(-pitch_KF_right_thigh,'minpeakheight',7,'minpeakdistance',200);

            [value_angle_pos, loc_angle_pos] = findpeaks(pitch_KF_right_thigh(...
                loc_angle_neg(1):length(pitch_KF_right_thigh)), 'minpeakheight',7, 'minpeakdistance',200,...
                'Npeaks',length(loc_angle_neg));

            loc_angle_pos = loc_angle_pos + loc_angle_neg(1);
            stride_GW (1,2) = mean(diff(loc_angle_pos));
            stride_GW (2,2) = abs( mean(diff(loc_angle_pos)-stride_GW (1,2)));

%             swing_GW (1,3) = mean(loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1));
%             swing_GW (2,3) =abs( mean((loc_angle_pos(1:length(loc_angle_pos)-1) ...
%                 - loc_angle_neg(1:length(loc_angle_pos)-1) )- swing_GW(1,3)));

            angle_GW(1,2) = mean(value_angle_pos);
            angle_GW(2,2) = mean(value_angle_neg);
            
            % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                input_signal,lwin_fsd,shift_fsd);
            
            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_right_thigh = gw.fusionGKF(gy,pitch_acc_right_thigh,...
                f,var(pitch_acc_right_thigh),var(pitch_acc_right_thigh),...
                var(gy),alpha1,alpha2,beta1,beta2,pitch_acc_right_thigh(1),...
                marker_fsd);
            
            if strcmpi(comparation,'yes')
                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_right_thigh(1);
                pitch_gyro_right_thigh = gw.integRate(1/f,gy,ini_pos);
                figure
                plot(time,pitch_acc_right_thigh)
                hold on
                plot(time,pitch_gyro_right_thigh,'r')
                plot(time,pitch_KF_right_thigh ,'black')
                plot(time,pitch_GKF_right_thigh,'cyan')
                title('Pitch angle of right thigh - Comparison')
                xlabel('Time (s)')
                ylabel('Pitch (deg)')
                legend('Accelerometer-based','Angular rate integration',...
                    'Kalman filter','Gated Kalman filter')   
%             else
%                 figure
%                 plot(time,pitch_KF_right_thigh ,'black')
%                 title('Pitch angle of right thigh')
%                 xlabel('Time (s)')
%                 ylabel('Pitch (deg)')
%                 legend('Kalman filter')
             end
            
        case 'left thigh'
            % Define the necessary signals.
            ax = a_X_left_thigh_1_C;
            az = a_Z_left_thigh_1_C;
            gy = g_Y_left_thigh_1_C;
            
            % Compute pitch using acceleration.
            pitch_acc_left_thigh = atan2d(az,ax);
            pitch_acc_left_thigh = gw.correct_quad_shifts(pitch_acc_left_thigh,'deg');
            pitch_acc_left_thigh = pitch_acc_left_thigh - 90;
            
            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_left_thigh = gw.fusion_KF(gy, pitch_acc_left_thigh,...
                f, var(pitch_acc_left_thigh) ,var(pitch_acc_left_thigh),...
                var(gy), alpha_KF, beta_KF, pitch_acc_left_thigh(1));
            
           % Centre the signal.
           pitch_KF_left_thigh = pitch_KF_left_thigh - pitch_KF_left_thigh(1);
           
           % Calculte features to characterise the movement.
            [value_angle_neg,loc_angle_neg] = findpeaks(-pitch_KF_left_thigh,'minpeakheight',7,'minpeakdistance',200);

            [value_angle_pos, loc_angle_pos] = findpeaks(pitch_KF_left_thigh(...
                loc_angle_neg(1):length(pitch_KF_left_thigh)), 'minpeakheight',0.5, 'minpeakdistance',200,...
                'Npeaks',length(loc_angle_neg));

            loc_angle_pos = loc_angle_pos + loc_angle_neg(1);
            stride_GW (1,4) = mean(diff(loc_angle_pos));
            stride_GW (2,4) = abs(mean(diff(loc_angle_pos)-stride_GW (1,4)));
% 
%             swing_GW (1,4) = mean(loc_angle_pos(1:length(loc_angle_pos)-1) - loc_angle_neg(1:length(loc_angle_pos)-1));
%             swing_GW (2,4) = abs( mean((loc_angle_pos(1:length(loc_angle_pos)-1)  - loc_angle_neg(1:length(loc_angle_pos)-1) )- swing_GW(1,4)));

            angle_GW(1,4) = mean(value_angle_pos);
            angle_GW(2,4) = mean(value_angle_neg);
           
            % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                input_signal,lwin_fsd,shift_fsd);
            
            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_left_thigh = gw.fusionGKF(gy,pitch_acc_left_thigh,...
                f,var(pitch_acc_left_thigh),var(pitch_acc_left_thigh),...
                var(gy),alpha1,alpha2,beta1,beta2,pitch_acc_left_thigh(1),...
                marker_fsd);
                       
            if strcmpi(comparation,'yes')
                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_left_thigh(1);
                pitch_gyro_left_thigh = gw.integRate(1/f,gy,ini_pos);
                figure
                plot(time,pitch_acc_left_thigh)
                hold on
                plot(time,pitch_gyro_left_thigh,'r')
                plot(time,pitch_KF_left_thigh ,'black')
                plot(time,pitch_GKF_left_thigh,'cyan')
                title('Pitch angle of left thigh - Comparison')
                xlabel('Time (s)')
                ylabel('Pitch (deg)')
                legend('Accelerometer-based','Angular rate integration',...
                    'Kalman filter','Gated Kalman filter')   
%             else
%                 figure
%                 plot(time,pitch_KF_left_thigh ,'black')
%                 title('Pitch angle of left thigh')
%                 xlabel('Time (s)')
%                 ylabel('Pitch (deg)')
%                 legend('Kalman filter')
            end
            
            
%         case 'left arm'
%             % Define the necessary signals.
%             gx = g_X_left_arm_1_C;
%             gy = g_Y_left_arm_1_C;
%             
%             % Integrate angular rate to compute pitch and roll.
%             pitch_gyro_left_arm = gw.integRate(1/f,gy,0);
%             roll_gyro_left_arm = gw.integRate(1/f,gx,0);
%             
%             % Apply high pass filter to partially remove bias.
%             lower_freq_limit=0.2; 
%             [b,a]=butter(3,lower_freq_limit*2/(f),'high'); 
%             pitch_gyro_left_arm_f = filtfilt(b,a,pitch_gyro_left_arm) + ...
%                 pitch_gyro_left_arm(1);
%             roll_gyro_left_arm_f  = filtfilt(b,a,roll_gyro_left_arm) + ...
%                 roll_gyro_left_arm(1);
            
%             % Plot results.
%             figure
%             plot(time,pitch_gyro_left_arm)
%             hold on
%             plot(time,pitch_gyro_left_arm_f,'black')
%             plot(time,roll_gyro_left_arm,'r')
%             plot(time,roll_gyro_left_arm_f,'g')
%             title('Pitch and roll angles of left arm')
%             legend('Pitch (unfiltered)', 'Pitch (filtered)',...
%                 'Roll (unfiltered)', 'Roll (filtered)');
%             
%         case 'right arm'
%             % Define the necessary signals.
%             gx = g_X_right_arm_1_C;
%             gy = g_Y_right_arm_1_C;
%             
%             % Integrate angular rate to compute pitch and roll.
%             pitch_gyro_right_arm = gw.integRate(1/f,gy,0);
%             roll_gyro_right_arm = gw.integRate(1/f,gx,0); 
%             
%             % Apply high pass filter to partially remove bias.
%             lower_freq_limit=0.2; 
%             [b,a]=butter(3,lower_freq_limit*2/(f),'high'); 
%             pitch_gyro_right_arm_f = filtfilt(b,a,pitch_gyro_right_arm)+...
%                 pitch_gyro_right_arm(1);
%             roll_gyro_right_arm_f  = filtfilt(b,a,roll_gyro_right_arm)+...
%                 roll_gyro_right_arm(1);
%             
%             % Plot results.
%             figure
%             plot(time,pitch_gyro_right_arm)
%             hold on
%             plot(time,pitch_gyro_right_arm_f,'black')
%             plot(time,roll_gyro_right_arm,'r')
%             plot(time,roll_gyro_right_arm_f,'g')
%             title('Pitch and roll angles of right arm')
%             legend('Pitch (unfiltered)', 'Pitch (filtered)',...
%                 'Roll (unfiltered)', 'Roll (filtered)');
%             
%         case 'center trunk'
%             % Define the necessary signals.
%             delay = 1;
%             ax = a_X_center_trunk_3_C(delay:end);
%             ay = a_Y_center_trunk_3_C(delay:end);
%             az = a_Z_center_trunk_3_C(delay:end);
%             gx = g_X_center_trunk_1_C(delay:end);
%             gy = g_Y_center_trunk_1_C(delay:end);
%             gz = g_Z_center_trunk_1_C(delay:end);
%             hx = h_X_center_trunk_3_C(delay:end);
%             hy = h_Y_center_trunk_3_C(delay:end);
%             hz = h_Z_center_trunk_3_C(delay:end);  
%             
%             % Compute roll, pitch using accelerometer.
%             pitch_acc = atan2d(sqrt(ay.^2+az.^2),ax) - 90;
%             roll_acc = atan2d(ay,az);
%                 
%             % Find the projections of hx and hy in the XY plane.
%             hx_n = -hx.*cosd(pitch_acc) + hy.*sind(pitch_acc).*...
%                 sind(roll_acc) - hz.*sind(pitch_acc).*cosd(roll_acc);
%             hy_n = hy.*cosd(roll_acc) + hz.*sind(roll_acc);
%                        
%             % Compute yaw.
%             yaw_mag = atan2d(hx_n,hy_n);
%             
%             % Compensate quadrante discontinuities.
%             yaw_mag = gw.correct_yaw_quad_shifts(yaw_mag,'deg');
%          
%             % Compute roll and pitch using non-quaternion Kalman Filter. 
%             alpha_KF = 100;
%             beta_KF = 0.01;
%             pitch_KF_center_trunk = gw.fusion_KF(gy,pitch_acc,f,...
%                 var(pitch_acc),var(pitch_acc),var(gy),alpha_KF,beta_KF,...
%                 pitch_acc(1));
%             roll_KF_center_trunk = gw.fusion_KF(gx,roll_acc,f,...
%                 var(roll_acc),var(roll_acc),var(gx),alpha_KF,beta_KF,...
%                 roll_acc(1));
%             alpha_KF = 1000;
%             beta_KF = 0.001;
%             yaw_KF_center_trunk = gw.fusion_KF(gz,yaw_mag,f,...
%                 var(yaw_mag),var(yaw_mag),var(gz),alpha_KF,beta_KF,...
%                 yaw_mag(1));
%             
%                         
%             % Find the projections of hx and hy in the XY plane using fused
%             % pitch and roll.
%             hx_n_KF = -hx.*cosd(pitch_KF_center_trunk') + hy.*sind(pitch_KF_center_trunk').*...
%                 sind(roll_KF_center_trunk') - hz.*sind(pitch_KF_center_trunk').*cosd(roll_KF_center_trunk');
%             hy_n_KF = hy.*cosd(roll_KF_center_trunk') + hz.*sind(roll_KF_center_trunk');
%             
%             % Compute yaw.
%             yaw_mag_kf = atan2d(hx_n_KF,hy_n_KF);
%             
%             % Compensate quadrante discontinuities.
%              yaw_mag_kf = gw.correct_yaw_quad_shifts(yaw_mag_kf,'deg');
%             correct_again = 1;
%             limit = 250;
%             while correct_again == 1
%                 correct_again = 0;
%                 for l = 2: length(yaw_mag_kf)
%                     if abs(yaw_mag_kf(l)-yaw_mag_kf(l-1)) > limit
%                         correct_again = 1;
%                     end
%                 end
%                 if correct_again == 1
%                     yaw_mag_kf = gw.correct_yaw_quad_shifts(yaw_mag_kf,'deg');
%                 end
%             end          
%             yaw_2KF_center_trunk = gw.fusion_KF(gz,yaw_mag_kf,f,...
%                 var(yaw_mag_kf),var(yaw_mag_kf),var(gz),alpha_KF,beta_KF,...
%                 yaw_mag_kf(1));
%             
%             % Compute intensity level.         
%             lwin_fsd = 20;    
%             threshold_fsd = 3;    
%             shift_fsd = 19;    
%             lambda = 30;
%             input_signal = sqrt(ax.^2+ay.^2+az.^2)';
%             [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
%                 threshold_fsd);
%             [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
%                 input_signal,lwin_fsd,shift_fsd);
%             
%             % Estimate Euler angles using Gated Kalman filter.
%             alpha1 = 100;
%             alpha2 = 10000;
%             beta1 = 0.001;
%             beta2 = 0.00001;
%             pitch_GKF_center_trunk = gw.fusionGKF(gy,pitch_acc,f,...
%                 var(pitch_acc),var(pitch_acc),var(gy),alpha1,alpha2,...
%                 beta1,beta2,pitch_acc(1),marker_fsd);
%             roll_GKF_center_trunk = gw.fusionGKF(gx,roll_acc,f,...
%                 var(roll_acc),var(roll_acc),var(gx),alpha1,alpha2,beta1,...
%                 beta2,roll_acc(1),marker_fsd);
%             yaw_GKF_center_trunk = gw.fusionGKF(gz,yaw_mag,f,var(yaw_mag),...
%                 var(yaw_mag),var(gz),alpha1,alpha2,beta1,beta2,yaw_mag(1),...
%                 marker_fsd);
%             yaw_2GKF_center_trunk = gw.fusionGKF(gz,yaw_mag_kf,f,var(yaw_mag_kf),...
%                 var(yaw_mag_kf),var(gz),alpha1,alpha2,beta1,beta2,yaw_mag_kf(1),...
%                 marker_fsd);
%             
%             % Compute roll, pitch and yaw using quaternion Extended Kalman 
%             % Filter.
%             gyroVarX = 0.1;
%             gyroVarY = 0.1;
%             gyroVarZ = 0.1;
%             mu_gain = 5;
%             alpha = 5;
%             q_ini = gw.eulerToQuat(roll_acc(1)/180*pi,pitch_acc(1)/180*pi,...
%                 yaw_mag(1)/180*pi);
% 
%             state_ini = [q_ini(1), q_ini(2), q_ini(3), 0.5]';
% 
%             [roll_EKF_center_trunk, pitch_EKF_center_trunk, ...
%                 yaw_EKF_center_trunk] = gw.quat9dofEKF(ax, ay, az, gx, gy,...
%                 gz, hx, hy, -hz, gyroVarX, gyroVarY, gyroVarZ, alpha,...
%                 mu_gain, f, state_ini);           
%              
%             if strcmpi(comparation,'yes')
% 
%                 pitch_gyro = gw.integRate(1/f,gy,pitch_acc(1));
%                 roll_gyro = gw.integRate(1/f,gx,roll_acc(1)); 
%                 yaw_gyro = gw.integRate(1/f,gz,yaw_mag(1)); 
%                 
%                 figure
%                 subplot(3,1,1)
%                 plot(time,pitch_acc)
%                 title('Orientation angles of trunk - Comparison')
%                 hold on
%                 plot(time,pitch_gyro,'r')
%                 plot(time,pitch_KF_center_trunk,'black')
%                 plot(time,pitch_GKF_center_trunk,'cyan')
%                 % plot(time,180/pi*pitch_EKF_center_trunk,'g')
%                 xlabel('Time (s)')
%                 ylabel('Pitch angle (deg)')
%                 legend('Accelerometer+Magnetometer','Gyroscope',...
%                     'Kalman filter','Gated Kalmal filter')
%                 subplot(3,1,2)
%                 plot(time,roll_acc)
%                 hold on
%                 plot(time,roll_gyro,'r')
%                 plot(time,roll_KF_center_trunk,'black')
%                 plot(time,roll_GKF_center_trunk,'cyan')
%                 %plot(time,180/pi*roll_EKF_center_trunk,'g')
%                 xlabel('Time (s)')
%                 ylabel('Roll angle (deg)')
%                 legend('Accelerometer+Magnetometer','Gyroscope',...
%                     'Kalman filter','Gated Kalman filter')
%                 subplot(3,1,3)
%                 plot(time,yaw_mag_kf)
%                 hold on
%                 plot(time,yaw_gyro,'r')
%                 plot(time,yaw_KF_center_trunk,'black')
%                 plot(time,yaw_GKF_center_trunk,'cyan')
%                 plot(time,yaw_2KF_center_trunk,'g')
%                 plot(time,yaw_2GKF_center_trunk,'--black')
%                 %plot(time,180/pi*yaw_EKF_center_trunk,'g')
%                 xlabel('Time (s)')
%                 ylabel('Yaw angle (deg)')
%                 legend('Accelerometer+Magnetometer','Gyroscope',...
%                     'Kalman filter','Gated Kalman filter',...
%                     'Double Kalman Filter','Double Gated Kalman Filter')      
%             else
%                 figure
%                 subplot(3,1,1)
%                 title('Orientation angles of trunk')
%                 plot(time,roll_KF_center_trunk)
%                 xlabel('Time (s)')
%                 ylabel('Roll angle (deg)')
%                 subplot(3,1,2)
%                 plot(time,pitch_KF_center_trunk)
%                 xlabel('Time (s)')
%                 ylabel('Pitch angle (deg)')
%                 subplot(3,1,3)
%                 plot(time,yaw_KF_center_trunk)
%                 xlabel('Time (s)')
%                 ylabel('Yaw angle (deg)')
%             end
    end
end

% ---------------------------------------------------------------------
% 2.9) Clear variables.
% ---------------------------------------------------------------------

clearvars -except filename_GW pitch_KF_right_shank pitch_KF_left_shank ...
    pitch_KF_right_thigh pitch_KF_left_thigh time Selection showPlots ...
    stride_GW swing_GW angle_GW

%--------------------------------------------------------------------------
% 3) Calculate pitch with Qualisys System.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 3.1) Load data.
% -------------------------------------------------------------------------

% Prompt to select the data file which is to be analyzed. 
% [filename, filepath] = uigetfile('*.mat', ...
%     'Select the QS data file(.mat)', '../../../Treadmill experiments/QS data');

% Search the appropiate QS file. To do that, we check the excel with the
% correspondences, so we can find the QS file from the GW file selected
% above.

% Read data from *.xlsx where are stored all filenames.
[~,file_excel] = xlsread( 'RecordingsGWandQSdata');
[rows,columns] = size(file_excel);

for i = 4:rows
    % Check if the cell is empty.
    if  isempty( char(strtrim(file_excel(i,9))) ) == 0  
        
        % We extract the GW name and find in the excel file.
        filename = strtrim(file_excel(i,9));

        if(filename{1,1} == filename_GW)
            filename_QS = strtrim(file_excel(i,6));
        end   
    end
end

% Load data.
load(fullfile('../../../Treadmill experiments/QS data/last recordings treadmill/',char( filename_QS)));

% Remove the '.mat' extension from the file name.
filename_QS = char(filename_QS);
filename = filename_QS(1:end-4);

% -------------------------------------------------------------------------
% 3.2) Read data.
% -------------------------------------------------------------------------

% Read sampling rate. The sampling rate is stored inside the 'FrameRate'
% field of structure which name is equal to the value of 'filename'.
SampRate = eval([filename, '.FrameRate']);

% Read markers data. Data are stored in the 'Trajectories.Labeled.Data'
% field.
data = eval([filename, '.Trajectories.Labeled.Data']);

% Remove the marker error information which is also stored by the Qualisys
% software. By doing this we only leave the (X, Y, Z) data.
data(:,4,:) = []; 

% Change the polarity of the data to match the desired coordinate frame.
data = data * -1;    

% Get size of the data matrix.
[n_markers, dummy, n_frames] = size(data);

% Extract data labels.
data_labels = eval([filename,'.Trajectories.Labeled.Labels']);

% Apply some slight modifications to the name of the labels. More
% especifically the underscores are removed and the substrings 'lo' and
% 'up' are replaced by 'lower' and 'upper' respectively.
for n = 1:n_markers
    data_labels{n} = strrep(data_labels{n},'_',' ');
    if ~isempty(strfind(data_labels{n},'lo'))
       data_labels{n} = strrep(data_labels{n},'lo','lower'); 
    end
    if ~isempty(strfind(data_labels{n},'up'))
       data_labels{n} = strrep(data_labels{n},'up','upper'); 
    end
end     

% Remove data structure for memory reasons.
eval(['clear ' filename]); 

% -------------------------------------------------------------------------
% 3.3) Interpolate data to remove NaN values (if existing).
% -------------------------------------------------------------------------
% The data matrix may have some 'NaN' values which are written by the
% Qualisys software in those instants in which the markers are not visible
% by any of the eight cameras. Markers are mainly ocluded by loose clothes
% so it is necessary to deal with this situation. We will apply
% interpolation to solve this issue.

% We first show how many data are missing for each one of the markers.
% nan_list = cell(1,n_markers);
% for n = 1:n_markers
%     nan_list{n} = [num2str(n),')  ', 'Marker: ',data_labels{n},'. Number of Nan values:', ...
%         num2str(sum(isnan(squeeze(data(n,1,:)))))];
% end
% 
% [Selection_correct,ok] = listdlg('ListString',nan_list,'Name',...
%     'Select the marker data you wish to correct','ListSize',[450 400],'SelectionMode',...
%     'multiple');

% calculate angular position of each marker trialgle in 2D
% and put them into data array Q_leg_pitch
% [pitch_KF_right_shank, pitch_KF_right_thigh, pitch_KF_left_shank, pitch_KF_left_thigh];

Q_leg_pitch = zeros(length(data),length(Selection));
time_vector = 1/200:1/200:length(data)/200;
if(length(time_vector)>length(time)) time_vector=time;end

order = {   'right lower shank','right upper shank','right back shank';...
            'right lower thigh','right upper thigh','right back thigh';...
            'left lower shank','left upper shank','left back shank';...
            'left lower thigh','left upper thigh','left back thigh';...
            'right hip','back hip','left hip'};

% Calculate the angle in the XZ plane.        
x=1;z=3;

for n=1:length(Selection)  % 4 leg segments (right shank, right thigh, left shank, left thigh)
    m1=0;m2=0;m3=0;
    
    for n2=1:3  % 3 markers (upper, lower and back)
        
        % Find the index for data vector.
        ind = find(strcmp(order{n,n2},data_labels));
        if isempty(ind)
            fprintf('Error: cannot find match for %s\n',order{n,n2})
            beep
        else
            switch n2
                case 1, m1 = ind;
                case 2, m2 = ind;
                case 3, m3 = ind;
            end
        end
    end
    
    
    
    seg_1 = -atan(squeeze((data(m2,x,:)-data(m1,x,:))./(data(m1,z,:)-data(m2,z,:))));        % marker 1,2: x/z
    seg_2 = -atan(squeeze((data(m2,z,:)-data(m3,z,:))./(data(m2,x,:)-data(m3,x,:))));        % marker 2,3: z/x
    seg_3 = -atan(squeeze((data(m1,z,:)-data(m3,z,:))./(data(m1,x,:)-data(m3,x,:))));        % marker 3,1: x/z
    
    % Calculate the angle in dregree.
     Q_leg_pitch(:,n) = (mean([seg_1 seg_2 seg_3],2).*180)./pi;
    
    % Centre the signals.
    Q_leg_pitch(:,n) = Q_leg_pitch(:,n) - Q_leg_pitch(1,n);
    
    % Calculte features to characterise the movement.
    [value_angle_neg,loc_angle_neg] = findpeaks(-Q_leg_pitch(:,n),'minpeakheight',7,'minpeakdistance',200);
    
    [value_angle_pos, loc_angle_pos] = findpeaks(Q_leg_pitch(...
        loc_angle_neg(1):length(Q_leg_pitch),n), 'minpeakheight',7, 'minpeakdistance',200,...
        'Npeaks',length(loc_angle_neg));
    
    loc_angle_pos = loc_angle_pos + loc_angle_neg(1);
    stride_QS (1,n) = mean(diff(loc_angle_pos));
    stride_QS (2,n) = abs( mean(diff(loc_angle_pos)-stride_QS (1,n)));
    
%     swing_QS (1,n) = mean(loc_angle_pos(1:length(loc_angle_pos)-1) - loc_angle_neg(1:length(loc_angle_pos)-1));
%     swing_QS (2,n) = abs(mean((loc_angle_pos(1:length(loc_angle_pos)-1)  - loc_angle_neg(1:length(loc_angle_pos)-1) )- swing_QS(1,n)));
    
    angle_QS(1,n) = mean(value_angle_pos);
    angle_QS(2,n) = mean(value_angle_neg);
    
end

% -------------------------------------------------------------------------
% 4) Comparation.
% -------------------------------------------------------------------------
% Others parameters.
mean_stride_GW_QS = mean([stride_QS(1,:) stride_GW(1,:)]);
mean_var_stride_GW = mean(stride_GW(2,:));
mean_var_stride_QS = mean(stride_QS(2,:));
mean_diff_angle = mean(abs([angle_GW(1,:) angle_GW(2,:)] -[angle_GW(1,:) angle_GW(2,:)]));
mean_diff_stride = mean(stride_QS(1,:) - stride_GW(1,:));

if strcmpi(showPlots,'yes')
for i=1:length(Selection)
    switch i
        case 1, 
            figure()
            plot(time_vector,Q_leg_pitch(1:length(time_vector),1));
            hold on
            plot(time_vector, pitch_KF_right_shank(1:length(time_vector)),'r');
            title('Pitch angle of right shank with Kalman filter')
            xlabel('Time (s)')
            ylabel('Pitch (deg)')
            legend('Qualisys System','Gait Watch')
            
        case 2,
            figure()
            plot(time_vector,Q_leg_pitch(1:length(time_vector),2));
            hold on
            plot(time_vector, pitch_KF_right_thigh(1:length(time_vector)),'r');
            title('Pitch angle of right thigh with Kalman filter')
            xlabel('Time (s)')
            ylabel('Pitch (deg)')
            legend('Qualisys System','Gait Watch')
            
        case 3, 
            figure()
            plot(time_vector,Q_leg_pitch(1:length(time_vector),3));
            hold on
            plot(time_vector, pitch_KF_left_shank(1:length(time_vector)),'r');
            title('Pitch angle of left shank with Kalman filter')
            xlabel('Time (s)')
            ylabel('Pitch (deg)')
            legend('Qualisys System','Gait Watch')
            
        case 4, 
            figure()
            plot(time_vector,Q_leg_pitch(1:length(time_vector),4));
            hold on
            plot(time_vector, pitch_KF_left_thigh(1:length(time_vector)),'r');
            title('Pitch angle of left thigh with Kalman filter')
            xlabel('Time (s)')
            ylabel('Pitch (deg)')
            legend('Qualisys System','Gait Watch')
            
    end
end
end

% save(['../../data/APA Parameters/GWvsQS/' ...
%        filename],'stride_GW','stride_QS', 'angle_GW','angle_QS',...
%        'mean_stride_GW_QS', 'mean_var_stride_GW', 'mean_var_stride_QS', ...
%        'mean_diff_angle', 'mean_diff_stride' );
