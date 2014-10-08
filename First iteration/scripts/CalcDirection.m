function [yaw_complete, north ] = CalcDirection(compass_trace,xyz_acc,samp_rate)
% CalcDirectionKB
% yaw_complete is in arcus not degree
% calculate direction of patient with GaitWatch.
% uses pitch and roll to correct data
% xyz_acc must be calibrated (-1 : +1)!
% Mag data must not be calibrated!
% samp_rate is needed to filter pitch and roll
% for test files see: GW Protocoll March 2014.txt
% yaw_complete is cumulative sum of direction changes
% Rotation left: increase of digits
% rotation right: decrease of numbers
% north is the reference, i. 0 means patient is heading north
% CalibrationList.mat should be on the Matlab path


% -------------------------------------------------------------------------
% 0) Initial configurataion.
% -------------------------------------------------------------------------
% clear all
% close all
% Load accelerometer's calibration parameters (bias_array, k_array)
load CalibrationList.mat      

% -------------------------------------------------------------------------
% 1) Load and extract data.
% -------------------------------------------------------------------------
% load GaitWatch_1301_020214_1055.mat   % rotation around z
% load C:\A-DATA\GaitWatchData\GaitWatch_1434_190414_2255.mat   % Rotation w. pitch and roll
% load GaitWatch_1433_190414_2142.mat   % Rotation with pitch
% load GaitWatch_1435_200414_1251.mat   % slow rotation for calibration north
% load GaitWatch_1418_150414_2256.mat   % file contains pitch in positive 
                                        % direction 360 deg plus roll in pos
                                        % direction 360 deg
                                        
% Extract magnetometer data channel.
% compass_trace = data(:,23);

% Get data length.
len_data = length(compass_trace);

% Get sampling rate.
% samp_rate = double(FileHeader.SampFreq);

% Extract acceleration data channels (x-y-z)
% xyz_acc = double([data(:,22) data(:,21) data(:,20)]);          
 
% -------------------------------------------------------------------------
% 2) Calibrate and clean acceleration data.
% -------------------------------------------------------------------------
% Apply bias (data obtained by Alberto).
% xyz_acc = xyz_acc + repmat(bias_array(22:-1:20)',len_data,1);   
% 
% % Apply scale factor (data obtained by Alberto).
% xyz_acc = xyz_acc .* repmat(k_array(22:-1:20)',len_data,1);    
%  
% Apply low pass filter.
upper_freq_limit = 1;
[b,a] = butter(3, upper_freq_limit * 2 / samp_rate);
xyz_acc = filtfilt(b,a,xyz_acc);

% -------------------------------------------------------------------------
% 3) Compute pitch and roll (orientation angles).
% -------------------------------------------------------------------------
roll   = atan2( xyz_acc(:,2) , sqrt(xyz_acc(:,1) .^2 + xyz_acc(:,3) .^2) );  
pitch  = atan2( -xyz_acc(:,1) , sqrt(xyz_acc(:,2) .^2 + xyz_acc(:,3) .^2) );
% pitch = filtfilt(b,a,pitch);
% this pitch is from -90 deg to +90 deg, it newer shows whwtehr the device
% is upside-down or not (180 deg)
% -------------------------------------------------------------------------
% 4) Process Magnetometer data.
% -------------------------------------------------------------------------
% Find Magnetometer-data
mag = cell(3,1);                                 

% Define channel's offsets.
offset1 = 20000;
offset2 = 32764;

% Select data using the offsets as thresholds.
mag{1} = find(compass_trace < -offset1/2);                                          
mag{2} = find(compass_trace >= -offset1/2 & compass_trace < offset1/2);
mag{3} = find(compass_trace >= offset1/2 & compass_trace < offset2);

% Get the length of the shortest array. 
min_val = min([length(mag{1}), length(mag{2}), length(mag{3})]);                    

% The follwing index is later used for displaying data vs time.
mag_time_ind = mag{1}(1:min_val);                                                   

% Parse data and crop data to the length of the shortest array. Offset is 
% not mag trace offset but serves to makk xyz data on one compass trace.
mag_data = double([compass_trace(mag{1}(1:min_val))+ offset1,...                    
                   compass_trace(mag{2}(1:min_val)),...
                   compass_trace(mag{3}(1:min_val))- offset1]);

% Remove spikes by interpolation.               
for n = 1:3                                                                                     
    mag_data(:,n) = spike_ex_2(mag_data(:,n),100);
end

% scale the mag traces

% scale = [830 970 630; -700 -780 -910];             % max and min of x-y-z compass traces (these are data I figured out, I am not sure they are ok)
% scale = [835 963 535; -657 -642 -905];               % here I use the min and max values which were obtained with a new online procedure.
% s2 = scale - repmat(mean(scale),2,1);                % scaling divisor
% bias_array(23:25) = - mean(scale);                   % new bias to be added to the raw values
% k_array (23:25) = 1.0 ./ s2(1,:);                    % new scaling factor
% save CalibrationList.mat bias_array k_array          % on this file we have Albertos scaling factors and bias data

% Calibrate magnetometer data
for n=1:3
    mag_data(:,n) = (mag_data(:,n) + bias_array(22+n)) * k_array (22+n);
end

% Rotate mag data with pitch and roll as correction (This is done to
% project the magnetometer data to the XY plane, where the atan can be
% applied to obtain heading.

% Last column stays empty. 
mag_data_rot = zeros(length(mag_data),4);             

for n=1:min_val   % Rotation around Y-axis.
    R = makehgtform('yrotate',pitch(mag_time_ind(n)),'xrotate',-roll(mag_time_ind(n)));
    mag_data_rot(n,:) = R * [mag_data(n,:),1]';
end

% Now the yaw (heading) can be computed using only X and Y data.
yaw = atan2(mag_data_rot(:,1),mag_data_rot(:,2));                                          

% Calculate sum of rotations to remove quadrant jumps.

% north =  1.79007906;
north = 2.01508879;
yaw(2:end,2) = diff(yaw);                                                          % yaw gets on more column 

ind = find((yaw(:,2) < -1) | (yaw(:,2) > 1));
for n=1:length(ind)
    for n2 = 1:10                                                               % search for last value that was not too big or small
        new = yaw(ind-n2,2);
        if abs(new) < 1
            yaw(ind,2) = new;
            break
        end
    end
end

yaw(1,2) = yaw(1,1);
yaw(:,2) = cumsum(yaw(:,2))-north;  

% Build time vector
% time =(0:len_data-1) / samp_rate;                 % time in s

% we would like to have the heading data at every data point, therefore:
% interpolation. Values at the begin and at the end are not precise if the
% patient moves in the begin and at the end of recording which should be
% omitted
% YY = spline(X,Y,XX)
yaw_complete = spline([1;mag_time_ind; len_data], [yaw(1,2);yaw(:,2);yaw(end,2)],1:len_data);



% -------------------------------------------------------------------------
% 5) Plot data.
% -------------------------------------------------------------------------
% Plot scaled mag data & pitch and roll
% figure
% subplot(2,1,1)
% plot(time(mag_time_ind),mag_data(:,1),'r')
% title('Mag Data')
% hold on
% plot(time(mag_time_ind),mag_data(:,2),'g')
% plot(time(mag_time_ind),mag_data(:,3),'b')
% axis tight
% legend('X axis', 'Y axis', 'Z axis')
% xlabel('Time (s)');
% ylabel('Magnetic field (Gauss)')
% 
% subplot(2,1,2)
% hold off
% plot(time,pitch,'r')
% hold on
% plot(time,roll,'g')
% axis tight
% legend('Pitch', 'Roll')
% xlabel('Time (s)');
% ylabel('Angle (rad)')
% 
% % Plot rotated and scaled mag data & pitch and roll.
% figure
% subplot(2,1,1)
% plot(time(mag_time_ind),mag_data_rot(:,1),'r')
% title('Rotated Mag Data')
% hold on
% plot(time(mag_time_ind),mag_data_rot(:,2),'g')
% plot(time(mag_time_ind),mag_data_rot(:,3),'b')
% axis tight
% 
% subplot(2,1,2)
% plot(time,pitch,'r')
% title('Pitch & Roll')
% hold on
% plot(time,roll,'g')
% axis tight
% 
% % Plot scaled mag data in 3D.
% figure
% plot3(mag_data(:,1), mag_data(:,2), mag_data(:,3)), axis equal
% xlabel('X'),ylabel('Y'),zlabel('Z')
% title('Scaled Mag data')
% 
% % Plot scaled and rotated mag data in 3D.
% figure
% plot3(mag_data_rot(:,1), mag_data_rot(:,2), mag_data_rot(:,3)), axis equal
% xlabel('X'),ylabel('Y'),zlabel('Z')
% title('Scaled & rotated Mag data')
% 
% % Plot yaw.
% figure
% subplot(2,1,1)
% plot(time(mag_time_ind),yaw(:,1));
% axis([time(1) time(end) -3.3 3.3])
% set(gca,'YTick',-pi:pi:pi,'YGrid','on','YLimMode','manual')
% title('yaw')
% ylabel('Radian')
% 
% subplot(2,1,2)
% ma= ceil(max(yaw(:,2))/(pi/2))*(pi/2);
% mi=floor(min(yaw(:,2))/(pi/2))*(pi/2);
% 
% plot(time(mag_time_ind),yaw(:,2));
% axis([time(1) time(end) mi ma])
% set(gca,'YTick',mi:pi/2:ma,'YGrid','on','YLimMode','manual','YTickLabel',{' '})
% line([time(1) time(end)],[0 0],'color',[1 0 0])
% ylabel('< r Rot [90 deg] l >')
% title('yaw cumulative')
% xlabel('Time [ s ]')
% 
% % plot(time(an:en),direction(an:en),'linewidth',1.5)
% % axis([axlim(1:2) mi ma])
% % set(gca,'YTick',mi:pi/2:ma,'YGrid','on','YLimMode','manual','YTickLabel',{' '})
% % % set(gca,'YTickLabel',{' '})
% % ylabel('[ 90 Deg ]')
% % title('Richtung')
end