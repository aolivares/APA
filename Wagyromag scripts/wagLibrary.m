function wag = wagLibrary

% FUNCTION ECLIBRARY contains the references to all necesarry functions to 
% load, extract, calibrate, process and plot data gathered using Wagyromag.
%
% -------------------------------------------------------------------------
% Author: Alberto Olivares.
% Entity: Universidad de Granada.
% Last modification: 20/03/2014.
% -------------------------------------------------------------------------

wag.load_data = @load_data;
wag.calibrate = @calibrate;
wag.pitch_roll_decomp = @pitch_roll_decomp;
wag.yaw_decomp = @yaw_decomp;

end
% END OF WAGLIBRARY FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [ax,ay,az,gx,gy,gz,hx,hy,hz,temp,bat,angle] = ...
    load_data(completepath)

% FUNCTION LOAD_DATA loads data from the .CSV files created by Wagyromag.
% The data are structured in the file as follows:
%   |_ Column 1: Acceleration (X axis).
%   |_ Column 2: Acceleration (Y axis).   
%   |_ Column 3: Acceleration (Z axis).
%   |_ Column 4: Angular rate (X axis).
%   |_ Column 5: Angular rate (Y axis).
%   |_ Column 6: Angular rate (Z axis).
%   |_ Column 7: Magnetic field (X axis).
%   |_ Column 8: Magnetic field (Y axis).
%   |_ Column 9: Magnetic field (Z axis).
%   |_ Column 10: Temperature.
%   |_ Column 11: Battery level.
%   |_ Column 12: Angle reference external input.
%
% *************************************************************************
% - Author: Alberto Olivares-Vicente. 
% - Institution: University of Granada, Spain.
% - Last revision: 20/03/2014.
% *************************************************************************
%
% Inputs:
%   |_'pathname' (string): Path to the folder where the data file is
%     stored.
%   |_'filename' (string): Name of the data file.
% 
% Outputs:
%   |_'ax' (vector):    Raw acceleration (X axis).
%   |_'ay' (vector):    Raw acceleration (Y axis).
%   |_'az' (vector):    Raw acceleration (Z axis).
%   |_'gx' (vector):    Raw angular rate (X axis).
%   |_'gy' (vector):    Raw angular rate (Y axis).
%   |_'gz' (vector):    Raw angular rate (Z axis).
%   |_'hx' (vector):    Raw magnetic field (X axis).
%   |_'hy' (vector):    Raw magnetic field (Y axis).
%   |_'hz' (vector):    Raw magnetic field (Z axis).
%   |_'temp' (vector):  Raw temperature.
%   |_'bat' (vector):   Raw battery level.
%   |_'angle' (vector): Raw reference angle.
% -------------------------------------------------------------------------

% Read the complete data matrix (all columns)
alldata = dlmread(completepath, ';', 3, 1); 

% Extract data columns.
ax = alldata(:,1);
ay = alldata(:,2);
az = alldata(:,3);
gx = alldata(:,4);
gy = alldata(:,5);
gz = alldata(:,6);
hx = alldata(:,7);
hy = alldata(:,8);
hz = alldata(:,9);
temp = alldata(:,10);
bat = alldata(:,11);
angle = alldata(:,12);

end
% END OF LOAD_DATA FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xC, yC, zC] = calibrate(x, y, z, id, sensor)

% FUNCTION CALIBRATE calibrates the raw data gathered by all three sensors
% (accelerometer, magnetometer and gyroscope).
% 
% *************************************************************************
% - Author: Alberto Olivares Vicente. 
% - Institution: University of Granada, Spain.
% - Last revision: 20/03/2014.
% *************************************************************************
%
% Inputs:
%   |_'x' (vector): Raw data along X-axis.
%   |_'y' (vector): Raw data along Y-axis.
%   |_'z' (vector): Raw data along Z-axis.
%   |_'id' (string): ID of the wagyromag.
%   |_'sensor' (string): Kind of sensor.
%
% Outputs:
%   |_'xC' (vector): Calibrated acceleration along X-axis.
%   |_'yC' (vector): Calibrated acceleration along Y-axis.
%   |_'zC' (vector): Calibrated acceleration along Z-axis. 
% -------------------------------------------------------------------------

% Check input parameters.
len = length(x);
if length(y)~=len || length(z)~=len
    error(['Input parameters ''X'', ''Y'' and ''Z'' must be ' ...
        'the same length.']);
end

if ~strcmpi(sensor,'acc') && ~strcmpi(sensor,'gyro') && ...
        ~strcmpi(sensor,'mag')
    error('The sensor name must be either ''acc'', ''gyro'' or ''mag''');
end

% Load Acceleromter, Magnetometer and Gyroscope calibration parameters.
load(strcat('data/calibration parameters/calib_params_wag',num2str(id),...
    '.mat'));

% Apply correction parameters.
xC = zeros(1, len);
yC = zeros(1, len);
zC = zeros(1, len);

% Selection of parameters according to the sensor.
if strcmpi(sensor,'acc')
    R = Ra;
    K = Ka;
    b = ba;
end
if strcmpi(sensor,'gyro')
    R = Rg;
    K = Kg;
    b = bg;
end
if strcmpi(sensor,'mag')
    R = Rh;
    K = Kh;
    b = bh;
end

for i = 1:len
    V = (inv(R) * inv(K) * ([x(i); y(i); z(i)] - b'))';
    xC(i) = V(1);  
    yC(i) = V(2);  
    zC(i) = V(3);
end

end
% END OF CALIBRATE FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


function [pitch, roll] = pitch_roll_decomp(ax, ay, az, units)

% FUNCTION PITCH_ROLL_ACC Computes pitch and roll angles using the
% decomposition of the acceleration measured with a triaxial accelerometer. 
% 
% *************************************************************************
% - Author: Alberto Olivares Vicente.
% - Entity: University of Granada, Spain. 
% - Last revision: 08/21/2013.
% *************************************************************************
%
% - Inputs:
%    |_'ax': vector containing the linear accelerations along the cartesian
%          X-axis.
%    |_'ay': vector containing the linear accelerations along the cartesian
%          Y-axis.
%    |_'az': vector containing the linear accelerations along the cartesian
%          Z-axis.
%    |_'units': string containing the units of angles to return. Set 'deg' 
%               or 'rad'.
%
% - Outputs:
%    |_pitch: vector containing the pitch calculated angle values given the
%             ay and az components of linear accelerations. Pitch is 
%             assumed to be the rotation angle about the Y-axis.
%    |_roll:  vector containing the roll calculated angle values given the 
%             ax and az components of linear accelerations. Roll is asumed 
%             to be the rotation angle about the X-axis.
% -------------------------------------------------------------------------

% Check input parameters.
if (length(ax)~=length(ay) || length(ay)~=length(az) || ...
    length(az)~=length(ax))
    error('Input vectors ax, ay and az must be the same length.');
end

if ~(strcmp(units, 'deg') || strcmp(units, 'rad'))
    error('Specified units are not correct. Set ''deg'' or ''rad''.');
end

% Compute pitch and roll. Pitch and roll values are obtained (in rad). Both
% pitch and roll values have different expressions depending on the 
% quadrant they are placed.
pitch = atan2(sqrt(ay.^2 + az.^2), -ax);
roll = atan2(az, ay);

% Roll and pitch values are transformed in degrees if so specified.
if strcmp(units, 'deg')
    roll = (180 / pi) .* roll;
    pitch = (180 / pi) .* pitch;
end

end
% END OF PITCH_ROLL_ACC FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [yaw, magXh, magYh] = yaw_decomp(pitch, roll, magX, magY, magZ)

% FUNCTION YAW_DECOMP computes the yaw angle using pitch and roll angles 
% and the magnetic field gathered with a triaxial magnetometer. 
%
% *************************************************************************
% - Authors: Alberto Olivares-Vicente and Gonzalo Ruiz-García.
% - Entity: University of Granada, Spain. 
% - Last revision: 20/03/2014.
% *************************************************************************

% - Inputs: 
%    |_ pitch: vector containing the a priori known pitch angle values.
%    |_ roll: vector containing the a priori known roll angle values.
%    |_ magX: vector containing the a priori known magnetic field component
%             along X-axis.
%    |_ magY: vector containing the a priori known magnetic field component
%             along Y-axis.
%    |_ magZ: vector containing the a priori known magnetic field component
%             along Z-axis.
%    |_ plotGraphics: string containing 'yes' or 'no' to specify if 
%                     results are to be plotted.
%
% - Outputs:
%    |_ yaw: vector containing the calculated yaw angles, considering both
%            local magnetic field and pitch and roll values.
%    |_ magXh: vector containing the projections of magX over the XY 
%              cartesian plane. These values are used to compute YAW output
%              values.
%    |_ magYh: vector containing the projections of magY over the XY 
%              cartesian plane. These values are used to compute YAW output
%              values.
% -------------------------------------------------------------------------

% Check input parameters.
len = length(pitch);
if (length(roll)~=len || length(magX)~=len || length(magY)~=len || ...
    length(magZ)~=len)
    error('All input arguments must be the same length.');
end

% Compute Yaw Angle. Calculation of horizontal magnetic X coordinate 
% component (magXh) and horizontal magnetic Y coordinate component (magYh).
% This is done by derotating the measured magnetic field in axes X and Y by
% the pitch and roll angles of the body. That is, we find the projection of
% the magnetic field in axes X and Y in the XY plane. 
magXh = magX.*cos(pitch) + magY.*sin(roll).*sin(pitch) - ...
        magZ.*cos(roll).*sin(pitch);
magYh = magY.*cos(roll) + magZ.*sin(roll);
yaw = atan((-magYh)./magXh);

% Quadrant compensations. This part may need to be changed depending on the
% body reference frame of the MIMU. 
for i=1:1:length(yaw)
    if magXh(i)>0 && magYh(i)==0
        yaw(i) = pi/2;
    end
    if magXh(i)<0 && magYh(i)==0
        yaw(i) =-pi/2;
    end
    if magXh(i)<0 && magYh(i)>0
        yaw(i) = yaw(i) - pi;
    end
    if magXh(i)<0 && magYh(i)<0
        yaw(i) = yaw(i) + pi;
    end
    
    while yaw(i) > pi 
        yaw(i) = yaw(i) - 2*pi; 
    end
    while yaw(i) <-pi 
        yaw(i) = yaw(i) + 2*pi; 
    end
end

end

% END OF YAW_DECOMP FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function quat_est = fqa(ax, ay, az, hx, hy, hz, magRef)

% FUNCTION FQA Applies the Factorized Quaternion Algorithm (without the 
% mechanism to avoid singularities) to estimate the orientation quaternions
% using the acceleration and magnetic field measured with a triaxial 
% accelerometer and a triaxial magnetometer.
% 
% *************************************************************************
% - Authors: Alberto Olivares-Vicente and Gonzalo Ruiz-García.
% - Entity: University of Granada, Spain. 
% - Last revision: 20/03/2014.
% *************************************************************************
%
% - Inputs:
%    |_ ax: vector containing linear acceleration along X-axis.
%    |_ ay: vector containing linear acceleration along Y-axis.
%    |_ az: vector containing linear acceleration along Z-axis.
%    |_ hx: vector containing magnetic component along X-axis.
%    |_ hy: vector containing magnetic component along Y-axis.
%    |_ hz: vector containing magnetic component along Z-axis.
%    |_ magRef: 3x1 matriz containing local magnetic reference vector.
%
% IMPORTANT NOTE: input parameters 'ax', 'ay', 'az', 'hx', 'hy'
%   and 'hz' must be the same length. Otherwise, an error is returned.
%
% - Outputs:
%    |_ quat: an Nx4 matrix containing calculated unit quaternions, where N
%             is length of input vectors.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0) CHECK INPUT PARAMETERS.
% -------------------------------------------------------------------------
len = length(ax);
if (length(ay)~=len || length(az)~=len || length(hx)~=len || ...
    length(hy)~=len || length(hz)~=len)
    error(['Input parameters ''ax'', ''ay'', ''az'', ''hx'', ' ...
        '''hy'' and ''hz'' must be the same length.']);
end
if (size(magRef, 1)~=3 || size(magRef, 2)~=1)
    error('Input parameter ''magRef'' must be a 3x1 real matrix.');
end

% -------------------------------------------------------------------------
% 1) Set variables.
% -------------------------------------------------------------------------
axn = zeros(len, 1);
ayn = zeros(len, 1);
azn = zeros(len, 1);
q_x = zeros(len, 4);
q_y = zeros(len, 4);
q_z = zeros(len, 4);

% Normalize magnetic reference.
magRefN = [magRef(1) magRef(2)]/norm([magRef(1) magRef(2)], 2);
quat_est = zeros(len, 4);

% -------------------------------------------------------------------------
% 2) Compute orientation quaternions.
% -------------------------------------------------------------------------
for i=1:1:len
    
    % Normalize acceleration.
    axn(i) = ax(i)/sqrt(ax(i).^2 + ay(i).^2 + az(i).^2);
    ayn(i) = ay(i)/sqrt(ax(i).^2 + ay(i).^2 + az(i).^2);
    azn(i) = az(i)/sqrt(ax(i).^2 + ay(i).^2 + az(i).^2);
    
    % Compute pitch quaternion.
    sen_pitch = axn(i); 
    cos_pitch = sqrt(1-sen_pitch.^2);
    sen_pitch_div2 = sign(sen_pitch)*sqrt((1-cos_pitch)/2);
    cos_pitch_div2 = sqrt((1+cos_pitch)/2);
    q_y(i, :) = cos_pitch_div2 * [1 0 0 0] + sen_pitch_div2 * [0 0 1 0];
        
    % Compute roll quaternion.
    sen_roll = -ayn(i)/cos_pitch;
    cos_roll = -azn(i)/cos_pitch;
    if sen_roll == 0 && cos_roll == -1
        sen_roll_div2 = sqrt((1-cos_roll)/2);
        cos_roll_div2 = sqrt((1+cos_roll)/2);
    else
        sen_roll_div2 = sign(sen_roll)*sqrt((1-cos_roll)/2);
        cos_roll_div2 = sqrt((1+cos_roll)/2);
    end
    q_x(i, :) = cos_roll_div2*[1 0 0 0]+sen_roll_div2*[0 1 0 0];
        
    % Compute yaw quaternion. 
        
        % Compute magnetic quaternion rotated to Earth's reference frame
        mag_quat = [0 hx(i) hy(i) hz(i)];
        mag_quat_rot = quat_mul(quat_inv(q_y(i,:)), ...
            quat_inv(q_x(i,:)));
        mag_quat_rot = quat_mul(mag_quat_rot, mag_quat);
        mag_quat_rot = quat_mul(mag_quat_rot, q_x(i,:));
        mag_quat_rot = quat_mul(mag_quat_rot, q_y(i,:));
        % Calculations.
        magN = [mag_quat_rot(2) mag_quat_rot(3)]/...
            norm([mag_quat_rot(2) mag_quat_rot(3)], 2);
        hxn = magN(1); 
        hyn = magN(2);
        aux = [hxn hyn; -hyn hxn] * [magRefN(1); magRefN(2)];
        cos_yaw = aux(1);
        sen_yaw = aux(2);
        sen_yaw_div2 = sign(sen_yaw)*sqrt((1-cos_yaw)/2);
        cos_yaw_div2 = sqrt((1+cos_yaw)/2);
        q_z(i, :) = cos_yaw_div2 * + [1 0 0 0] + sen_yaw_div2 * [0 0 0 1];
        
    % Compute final estimated quaternion. 
    quat_est(i, :) = quat_mul(q_z(i,:), q_y(i,:));
    quat_est(i, :) = quat_mul(quat_est(i,:), q_x(i,:));
    % Normalize final estimated quaternion. 
    quat_est(i, :) = quat_est(i, :) / norm(quat_est(i, :), 2);
    
end

end