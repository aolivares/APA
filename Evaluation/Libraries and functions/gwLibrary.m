function gw = gwLibrary

% FUNCTION GWLIBRARY contains the references to all necesarry functions to 
% load, extract, calibrate, process and plot data gathered using GaitWatch.
%
% -------------------------------------------------------------------------
% Author: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 13/12/2013.
% -------------------------------------------------------------------------
%
%  ******************* DATA FILE HANDLING FUNCTIONS ***********************
gw.openGWfile = @openGWfile;
gw.getDataChannel = @getDataChannel;
gw.getFHinfo = @getFHinfo;

% **************** MAGNETOMETER PREPROCESSING FUNCTIONS *******************
gw.getMagData = @getMagData;
gw.correct_mag_data = @correct_mag_data;

% ********************** CALIBRATION FUNCTIONS ****************************
gw.mag3DCal = @mag3DCal;
gw.acc1DCal = @acc1DCal;
gw.acc3DCal = @acc3DCal;
gw.gyro1DCal = @gyro1DCal;

% ********** COMPUTATION OF CALIBRATION PARAMETERS FUNCTIONS **************
gw.gyro1DscaleF = @gyro1DscaleF;
gw.comp_acc_mag_cal_params = @comp_acc_mag_cal_params;

% ********************** AUXILIARY FUNCTIONS ******************************
gw.getDCindexes = @getDCindexes;
gw.selectStaticPositions3D = @selectStaticPositions3D;
gw.get_acc_parallel_values = @get_acc_parallel_values;
gw.quatTOeuler = @quatTOeuler;
gw.eulerToQuat = @eulerToQuat;
gw.quat_mul = @quat_mul;
gw.plotData = @plotData;

% ****************** ORIENTATION ESTIMATION FUNCTIONS *********************
gw.quat9dofEKF = @quat9dofEKF;
gw.optimize_KF = @optimize_KF;
gw.fusion_KF = @fusion_KF;
gw.optimize_EKF = @optimize_EKF;
gw.correct_quad_shifts = @correct_quad_shifts;
gw.correct_yaw_quad_shifts = @correct_yaw_quad_shifts;
gw.correct_ekf_quad_shifts = @correct_ekf_quad_shifts;
gw.integRate = @integRate;
gw.fusionGKF = @fusionGKF;
% ********************** MOTION INTENSITY DETECTION ***********************
gw.fsd = @fsd;
gw.compEstMark = @compEstMark;

% **************************** MISCELLANEA ********************************
gw.sync_GW_QS = @sync_GW_QS;

end
% END OF GWLIBRARY FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **                   DATA FILE HANDLING FUNCTIONS                      **
% *************************************************************************

function [data, FileHeader] = openGWfile(dataPath)

% FUNCTION OPENGWFILE loads GaitWatch data file from the hard drive.
%
% Input parameters:
% |_ 'calPath':  Path to the file containing the GW calibration parameters.
% |_ 'dataPath': (optional) User-defined path to the GW data file. When
%                this parameter is not specified, the function will prompt 
%                the user to select the location of the data file.
%
% Output parameters:
% |_ 'data':      Matrix in which measured magnitudes are stored in 
%                 columns.
% |_ 'time':      Vector containing the time samples.
% |_ 'FH':        File header containing information relative to the 
%                 monitoring session.
% 
% -------------------------------------------------------------------------
% Authors: Alberto Olivares and Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Load data (either prompting the user or from the predefined path).

% If the user provides the data path, then load data from such data path.
if nargin==1            
    load(dataPath);
    
% If, on the other hand, the function is called without any arguments, then
% a file browser is prompted to the user to select the file to be loaded.    
elseif nargin==0  
    if ispc
        [filename, filepath] = uigetfile('\data\*.mat',...
            'Select GaitWatch data file');
    else
        [filename, filepath] = uigetfile('/data/*.mat',...
            'Select GaitWatch data file');
    end
    load(fullfile(filepath,filename));
end

end
% END OF OPENGWFILE FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function channel = getDataChannel(data_struct,magnitude, axis,...
    position, segment)

% FUNCTION GETDATACHANNEL returns the number of the channel associated to 
% the segment provided by the user. 
%
% Input parameters:
% |_ 'data_struct': GaitWatch's data struct (cell).
% |_ 'magnitude': Physical magnitude measured by the sensor (char).
% |_ 'axis': Axis of the sensor (char).
% |_ 'position': Position in which the sensor is place (string). 
% |_ 'segment': Segment from which we want to find the data channel 
%    (string). 
%
% Output parameters:
% |_ 'channel': Number of the channel in which the desired data can be
% found inside the 'data' matrix. 
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares and Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% The following loop sweeps 'data_struct' and extracts the channel number
% for the corresponding magnitude, axis and position. 
for i = 1: length(data_struct)
    if strcmpi(data_struct(i,2),magnitude) && strcmpi(data_struct(i,3),...
            axis) && strcmpi(data_struct(i,4),position) && ...
            strcmpi(data_struct(i,5),segment)
        channel = data_struct{i,1};
    end
end

end
% END OF GETDATA CHANNEL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [f, date, start_time, end_time, file_id] = ...
    getFHinfo(FileHeader)

% FUNCTION GETFHINFO gets information from the data file header.
%
% Input parameters:
% |_ 'FileHeader':  Header of the GW data file. The header contains
%                   different useful information about the data.
%
% Output parameters:
% |_ 'f':           Sampling frequency (scalar). 
% |_ 'date':        Date in the format dd/mm/yy (string).
% |_ 'start_time':  Start time in the format hh:mm:ss (string).
% |_ 'end_time':   	End time in the format hh:mm:ss (string).
% |_ 'file_id':     Identification number of the data file (scalar).
% |_ 'FileName':    File name given by the user.
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Extract sampling frequency and transform it from 'int16' to 'double'. 
f = double(FileHeader.SampFreq);

% 2) Extract initial and final time of the data gathering maneuvers. 
time_begin = double(FileHeader.time_begin);
time_end = double(FileHeader.time_end);
date = strcat(num2str(time_begin(5)),'/',num2str(time_begin(6)),'/',... 
    num2str(time_begin(7)));
start_time = strcat(num2str(time_begin(3)),':',num2str(time_begin(2)),...
    ':',num2str(time_begin(1)));
end_time = strcat(num2str(time_end(3)),':',num2str(time_end(2)),':',... 
    num2str(time_end(1)));

% 3) Extract file ID and file name.
file_id = double(FileHeader.FileNumber);
%FileName = FileHeader.FileName;

end
% END OF GETFHINFO FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **               MAGNETOMETER PREPROCESSING FUNCTIONS                  **
% *************************************************************************

function [hx, hy, hz] = getMagData(magnetometer)

% FUNCTION GETMAGDATA Separates the triaxial magnetometer data contained in
% the same data channel. Data from each axis has a different artificially 
% added offset, this way they can be easily separated.
% 
% Input parameters: 
% |_ 'magnetometer': Magnetometer signal containing data from all three
%                    axes (vector).
% Output parameters:
% |_ 'hx': Magnetometer data gathered along X axis (vector).
% |_ 'hy': Magnetometer data gathered along Y axis (vector).
% |_ 'hz': Magnetometer data gathered along Z axis (vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 05/11/2013.
% -------------------------------------------------------------------------

% 1) Set artificial offsets (these offsets are configured in the firmware
%    of GaitWatch).
offset1 = 20000;
offset2 = 32764;

% 2) Separate data. 
hx = magnetometer(magnetometer < -offset1/2) + offset1;
hy = magnetometer(magnetometer >= -offset1/2 & magnetometer < offset1/2); 
hz = magnetometer(magnetometer >= offset1/2 & magnetometer < offset2) - ...
    offset1; 

% 3) Crop the three signals to the length of the shortest signal.
len_hx = length(hx);
len_hy = length(hy);
len_hz = length(hz);

min_val = min([len_hx len_hy len_hz]);

hx = double(hx(1:min_val));
hy = double(hy(1:min_val));
hz = double(hz(1:min_val));

end
% END OF FUNCTION GETMAGDATA

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function h_corr = correct_mag_data(h)

% FUNCTION CORRECT_MAG_DATA removes any erroneous value in the magnetometer
% signal. The user is shown a figure depicting the magnetometer signal and 
% he should select (using the datacursor which is activated by default) the 
% data points (if existing) that need to be removed. When ready, the user 
% should click in the "Continue execution" button. 
% 
% Input parameters:
% |_ 'h': Magnetometer signal (vector).
%
% Output parameters:
% |_ 'h_corr': Corrected magnetometer signal (vector).
% 
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Initialize the value of the corrected magnetometer signal.
h_corr = h;

% 2) Define the title of the figure which will be shown to the user.
fig_title = 'SELECT THE ERRONEOUS VALUES TO BE REMOVED';

% 3) The signal is shown to the user so he can select (using the 
%    'datacursor') the erroneous data points. The indexes of these points 
%    are stored in the 'indexes' vector.
indexes = getDCindexes(h,fig_title);

% 4) Substitute the erroenous values by the subsequent value (unless the
%    erroneous value is in the last position of the singal, in which case 
%    it is substituted by the preceding value). 
for i = 1:length(indexes)
    if indexes(i) == length(h)
        h_corr(indexes(i)) = h(indexes(i)-1);
    else
        h_corr(indexes(i)) = h(indexes(i)+1);
    end
end

end
% END OF CORRECT_MAG_DATA FUNCTION 

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **                      CALIBRATION FUNCTIONS                          **
% *************************************************************************

function [cal_hx, cal_hy, cal_hz] = mag3DCal(raw_hx, raw_hy, raw_hz,...
    segment)

% FUNCTION MAG3DCAL calibrates raw triaxial magnetometer data. To do so, it
% loads the already existing calibration parameters from the hard drive.
%
% Input parameters:
% |_ 'raw_hx':  Raw magnetic field measured along X axis (N x 1 vector).
% |_ 'raw_hy':  Raw magnetic field measured along Y axis (N x 1 vector).
% |_ 'raw_hz':  Raw magnetic field measured along Z axis (N x 1 vector).
% |_ 'segment': Segment to which the unit containing the magnetometer is
%               attached (string).
% 
% Output parameters:
% |_ 'cal_hx': Calibrated magnetic field measured along X axis 
%              (N x 1 vector).
% |_ 'cal_hy': Calibrated magnetic field measured along Y axis 
%              (N x 1 vector).
% |_ 'cal_hz': Calibrated magnetic field measured along Z axis 
%              (N x 1 vector).
% 
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Build name of the file containing the calibration parameters.
cal_file_name = strcat(segment,'_XYZ_magCalParams');

% 2) Build name of the path to the file containing the calibration 
%    parameters.
data_path = strcat('Calibration Parameters/magnetometer/',cal_file_name);
    
% 3) Load calibration parameters.
try
    load(data_path)
catch
    error(['ERROR: Either the file containing the magnetometer ',...
        'calibration parameters does not exist, or the path is not ',...
        'correct']);
end

% 4) Initialize calibrated signals.
cal_hx = zeros(1,length(raw_hx));
cal_hy = zeros(1,length(raw_hy));
cal_hz = zeros(1,length(raw_hz));

% 5) Apply calibration parameters.
for i = 1:length(raw_hx)
    H = inv(alpha_h) * ([raw_hx(i);raw_hy(i);raw_hz(i)] - beta_h');
    cal_hx(i) = H(1);
    cal_hy(i) = H(2);
    cal_hz(i) = H(3);
end

end
% END OF MAG3DCAL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function cal_acc = acc1DCal(raw_acc, sens_axis, position, segment)

% FUNCTION ACC1DCAL calibrates raw uni-axial acceleration data.
%
% Input parameters:
% |_ 'raw_acc':   Raw acceleration gathered along a single axis (N x 1
%                 vector).
% |_ 'sens_axis': Label of the axis in which the acceleration was gathered 
%                 (char).
% |_ 'position':  Position of the unit which includes the accelerometer
%                 (string). 
% |_ 'segment':   Segment to which the unit which includes the 
%                 accelerometer is attached. 
%
% Output parameters:
% |_ 'cal_acc': Calibrated acceleration (N x 1 vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Build name of the file containing the calibration parameters.
cal_file_name = strcat(position,'_',segment,'_',sens_axis,'_accCalParams');

% 2) Build name of the path to the file containing the calibration 
%    parameters.
data_path = strcat('Calibration Parameters/accelerometer/',cal_file_name);

% 3) Load file containing the calibration parameters.
load(data_path);

% 4) Apply calibration parameters.
cal_acc = (raw_acc * k) + bias;

end
% END OF ACC1DCAL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [ax_cal, ay_cal, az_cal] = acc3DCal(ax_raw, ay_raw, az_raw,...
    segment)

% FUNCTION ACC3DCAL calibrates raw triaxial acceleration.
%
% Input parameters:
% |_ 'ax_raw':  Raw acceleration measured along X axis (N x 1 vector).
% |_ 'ay_raw':  Raw acceleration measured along Y axis (N x 1 vector).
% |_ 'az_raw':  Raw acceleration measured along Z axis (N x 1 vector).
% |_ 'segment': Segment to which the unit containing the accelerometer is
%               attached (string). 
%
% Output parameters:
% |_ 'ax_cal': Calibrated acceleration along X axis (N x 1 vector). 
% |_ 'ay_cal': Calibrated acceleration along Y axis (N x 1 vector). 
% |_ 'az_cal': Calibrated acceleration along Z axis (N x 1 vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 11/11/2013.
% -------------------------------------------------------------------------

% 1) Build name of the file containing the calibration parameters.
cal_file_name = strcat(segment,'_XYZ_accCalParams');

% 2) Build name of the path to the file containing the calibration 
%    parameters.
data_path = strcat('Calibration Parameters/accelerometer/',cal_file_name);

% 3) Load file containing the calibration parameters.
load(data_path);

% 4) Apply calibration parameters.
ax_cal = zeros(1,length(ax_raw));
ay_cal = zeros(1,length(ay_raw));
az_cal = zeros(1,length(az_raw));
for i=1:length(ax_raw)
    A = inv(alpha_acc) * ([ax_raw(i);ay_raw(i);az_raw(i)] - beta_acc');
    ax_cal(i) = A(1);
    ay_cal(i) = A(2);
    az_cal(i) = A(3);
end

end
% END OF ACC3DCAL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function cal_rate = gyro1DCal(raw_rate, sens_axis, position, segment)

% FUNCTION ACC1DCAL calibrates raw uniaxial angular rate data.
%
% Input parameters:
% |_ 'raw_rate':  Raw angular rate gathered along a single axis (N x 1
%                 vector).
% |_ 'sens_axis': Label of the axis in which the angular rate was gathered 
%                 (char).
% |_ 'position':  Position of the unit which includes the gyroscope
%                 (string). 
% |_ 'segment':   Segment to which the unit which includes the 
%                 gyroscope is attached. 
%
% Output parameters:
%
% |_ 'cal_rate': Calibrated angular rate (N x 1 vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Build name of the file containing the calibration parameters.
cal_file_name = strcat(position,'_',segment,'_',sens_axis,'_gyroCalParams');

% 2) Build name of the path to the file containing the calibration 
%    parameters.
data_path = strcat('Calibration Parameters/gyroscope/',cal_file_name);

% 3) Load file containing the calibration parameters.
load(data_path);

% 4) Apply calibration parameters.
cal_rate = (raw_rate - bias) / k_av;

end
% END OF GYRO1DCAL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **           COMPUTATION OF CALIBRATION PARAMETERS FUNCTIONS           **
% *************************************************************************

function k = gyro1DscaleF(time,wg,rotation)

% FUNCTION GYRO1DSCALEF computes the scale factor of a gyroscope axis given
% the raw angular rate gathered during a known rotation. 

% Input parameters:
% |_ 'time':      Time vector associated to the gathered rotation (vector). 
% |_ 'wg':        Raw angular rate (vector).
% |_ 'rotation':  Known rotation angle (scalar). 
% 
% Output parameters:
% |_ 'k': Computed scale factor. 
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) The raw angular rate is integrated and then divided by the known 
%    rotation angle to find the scale factor. 
int_wg = 0;
for i = 2:length(wg)
    int_wg = int_wg + 0.5 * (time(i) - time(i-1)) * (wg(i) + wg(i-1));
end

k = int_wg / rotation;
end
% END OF GYRO1DSCALEF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [alpha, beta] = comp_acc_mag_cal_params(sensor, x_axis, y_axis,...
    z_axis, ref_val)

% FUNCTION COMP_ACC_MAG_CAL_PARAMS Computes the calibration parameters of a
% triaxial accelerometer or a triaxial magnetoemter using an ellipsoid 
% fitting algorithm. 
%
% Input parameters:
% |_ 'sensor': Kind of sensor that we wish to calibrate ('acc' or 'mag').
% |_ 'x_axis': Signal containing the measured data along X axis.
% |_ 'y_axis': Signal containing the measured data along Y axis.
% |_ 'z_axis': Signal containing the measured data along Z axis.
% |_ 'ref_val': 
%
% Output parameters:
% |_ 'alpha': Matrix containing the scale factors (in its diagonal) and the
%             non orthogonality factors (3x3 matrix).
% |_ 'beta':  Vector containing the bias of each of the three axes (1x3
%             vector).
% 
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 06/11/2013.
% -------------------------------------------------------------------------

% 1) We first need to define some global variables as they will be used in 
%    a different function which does not take any parameters.
global x;
global y;
global z;
global val;

% 2) We assign each of the data signals to the global variables as well as 
%    the reference value (which is the magnitude of the Earth's gravity 
%    vector if we are calibrating the accelerometer or the magnitude of the
%    local Earth's magnetic field if we are calibrating the magnetometer).
x = x_axis;
y = y_axis;
z = z_axis;
val = ref_val;

% 3) Set the initial values of the calibration parameters according to the
%    sensor being calibrated. 
if strcmp(sensor,'mag')
    p0=[1e-3,1e-3,1e-3,0,0,0,80,20,-20];
elseif strcmp(sensor,'acc')
    p0=[1/1.8,1/1.8,1/1.8,1/1.8,1/1.8,1/1.8,1,1,1];
end

% 4) We now call the optimization process which is based in an ellipsoid
%    fitting procedure. The function to be optimized is defined in the
%    'ellipsoidEOF' function which can also be found in this library. The
%    estimation is done using Matlab's 'lsqnonlin' function which applies 
%    the Gauss-Newton algorithm to minimize a nonlinear function. 
p = lsqnonlin(@ellipsoidEOF, p0, -Inf, Inf, optimset('maxf', 500000,...
    'TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 1e6));

% 5) When the estimation of the parameters is over, we proceed to extract the
%    values from the vector. 

alphax = 1 / p(1);  % Scale factor (X axis).
alphay = 1 / p(2);  % Scale factor (Y axis).
alphaz = 1 / p(3);  % Scale factor (Z axis).
sxy = p(4);         % Non-orthogonality factor between X and Y axes. 
sxz = p(5);         % Non-orthogonality factor between X and Z axes. 
syz = p(6);         % Non-orthogonality factor between Y and Z axes. 
betax = p(7);       % Bias (X axis).
betay = p(8);       % Bias (Y axis).
betaz = p(9);       % Bias (Z axis).

% 6) We now build the matrix containing the scale factors and the
%    non-orthogonality factors as well as the bias vector. 
alpha=[alphax,sxy,sxz; sxy,alphay,syz; sxz,syz,alphaz];
beta=[betax, betay, betaz];

end
% END OF COMP_ACC_MAG_CAL_PARAMS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = ellipsoidEOF(p)

% FUNCTION ELLIPSOIDEOF defines the error function to be minimized by the 
% Gauss-Newton procedure which is called from 'comp_acc_mag_cal_params' 
% function. 
%
% Input parameters:
% |_ 'p': vector containing the calibration parameters (1 x 9 vector).
% 
% Output parameters:
% |_ 'F': value of the error function (scalar). 
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Definition of the global variables which will be used in this 
%    function.
global x
global y
global z
global val

% 2) Extraction of the calibration parametes from the input parameter 
%    vector.
alphax = p(1);  % Inverse scale factor (X axis).
alphay = p(2);  % Inverse scale factor (Y axis).
alphaz = p(3);  % Inverse scale factor (Z axis).
sxy = p(4);     % Non-orthogonalization factor between axes X and Y.
sxz = p(5);     % Non-orthogonalization factor between axes X and Z.
syz = p(6);     % Non-orthogonalization factor between axes Y and Z.
betax = p(7);   % Bias (X axis).
betay = p(8);   % Bias (Y axis).
betaz = p(9);   % Bias (Z axis).

% 3) Check GaitWatch's user manual to learn more about the error function.
F = val^2-((x-betax).*((x-betax).*(alphax^2+sxy^2+sxz^2)+(y-betay).*...
    (sxy*alphax+alphay*sxy+syz*sxz)+(z-betaz).*(sxz*alphax+syz*sxy+...
    alphax*sxz))+(y-betay).*((x-betax).*(alphax*sxy+sxy*alphay+sxz*syz)+...
    (y-betay).*(sxy^2+alphay^2+syz^2)+(z-betaz).*(sxz*sxy+syz*alphay+...
    alphaz*syz))+(z-betaz).*((x-betax).*(alphax*sxz+sxy*syz+sxz*alphaz)+...
    (y-betay).*(sxy*sxz+alphay*syz+syz*alphax)+(z-betaz).*(sxz^2+syz^2+...
    alphaz^2)));
end
% END OF ELLIPSOIDEOF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **                        AUXILIARY FUNCTIONS                          **
% *************************************************************************
 
function indexes = getDCindexes(signal,fig_title)

% FUNCTION GETDCINDEXES extracts the indexes of the data points shown in a 
% figure and selected using the data cursor tool.
%
% Input parameters:
% |_ 'signal': The signal to be plotted.
% |_ 'fig_title': Title which is shown in the figure. 
%
% Output parameters:
% |_ 'indexes': The selected indexes.
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Create figure displaying the signal.
fig = figure; 
plot(signal)
set(gcf,'position',[40 80 1600 900])
title(fig_title);

% 2) The data cursor mode is activated. This means that the data cursor 
%    tool is already active when the figure is shown and that the selected
%    data points are stored. The execution of the routine does not 
%    continue until the user clicks the "Continue execution" button. 
dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip',...
'SnapToDataVertex','off','Enable','on')
uicontrol('Style', 'pushbutton','String', 'Continue execution',...
    'Position', [720 20 150 20],'Callback', 'uiresume(gcbf)');
uiwait(fig);

% 3) Extract the indexes of the datacursors. "c_info" contains a structure
%    for each one of the selected data points. The indexes are contained
%    in the "DataIndex" field of the structure. 
c_info = getCursorInfo(dcm_obj);
indexes = zeros(1,length(c_info));
for i = 1:length(c_info)
    indexes(i) = c_info(i).DataIndex;
end
indexes = sort(indexes);
end
% END OF GETDCINDEXES FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


function [detected, xo, yo, zo] = selectStaticPositions3D(x, y, z, lmin,...
    dmax, dmin, showPlot)

% FUNCTION SELECTSTATICPOSITIONS3D Extracts the mode of each one of a set 
% of detected quasi-static intervals. The detection of these intervals is 
% controlled by the input parameters. 
%
% Input parameters:
% |_ 'x':       Data gathered along X axis (vector).
% |_ 'y':       Data gathered along Y axis (vector).
% |_ 'z':       Data gathered along Z axis (vector).
% |_ 'lmin':    Minimum length of quasi-static period (scalar).
% |_ 'dmax':    Maximum permited deviation from the reference value in a
%               quasi-static period (scalar).
% |_ 'dmin':    Minimum distance between consecutive intervals so they are
%               considered to be different (scalar).
%
% Output parameters:
% |_ 'detected':    Number of detected quasi-static periods (scalar).
% |_ 'xo':          Mode of the 'x' signal during each quasi-static period
%                   (1 x N where N is the number of detected periods).
% |_ 'yo':          Mode of the 'y' signal during each quasi-static period
%                   (1 x N where N is the number of detected periods).
% |_ 'zo':          Mode of the 'z' signal during each quasi-static period
%                   (1 x N where N is the number of detected periods).
%
% -------------------------------------------------------------------------
% Authors: Rafael L?pez, Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Find the static intervals complying with the provided conditions.
intervals = zeros(2,0);
newInterval = 1;
xref = 0;
yref = 0;
zref = 0;
pStart = 0;
pFinal = 0;
for i = 1 : length(x)
    if newInterval == 1
        xref = x(i);
        yref = y(i);
        zref = z(i);
        pStart = i;
        pFinal = pStart;
        newInterval = 0;
    else
        if abs(xref - x(i)) < dmax && abs(yref - y(i)) < dmax && ...
                abs(zref - z(i)) < dmax
            pFinal = pFinal + 1;
             if i == length(x) && (pFinal - pStart + 1) > lmin 
                intervals = [intervals [pStart; pFinal]];
                newInterval = 1;
                i = i - 1;
             end
        else
            if (pFinal - pStart + 1) < lmin
                i = pStart + 1;
            else
                intervals = [intervals [pStart; pFinal]];
                i = i - 1;
            end
            newInterval = 1;
        end
    end
end
interval_size = size(intervals);
detected = interval_size(2);

% 2) Once the intervals have beeen found, we compute the mode of each one 
%    of them for the three signals. 
if detected > 0
    xom = zeros(length(intervals(1,:)),1);
    yom = zeros(length(intervals(1,:)),1);
    zom = zeros(length(intervals(1,:)),1);

for i=1:length(intervals(1,:))
    xom(i) = mode(x(intervals(1,i) : intervals(2,i)));
    yom(i) = mode(y(intervals(1,i) : intervals(2,i)));
    zom(i) = mode(z(intervals(1,i) : intervals(2,i)));
end

interval = zeros(2,0);
for i = 1: length(intervals(1,:))
    if i == 1
        xo = xom(1);
        yo = yom(1);
        zo = zom(1);
        interval = intervals(:,1);
        n = 1;
    else
        if abs(xo(n) - xom(i)) > dmin || abs(yo(n) - yom(i)) > dmin ||...
                abs(zo(n) - zom(i)) > dmin
            xo = [xo; xom(i)];
            yo = [yo; yom(i)];
            zo = [zo; zom(i)];
            interval = [interval, intervals(:,i)];
            n = n + 1;
        end
    end
end
detected = length(xo);

% 3) Plot results (signals vs. detected quasi-static periods).
if strcmpi(showPlot,'yes')
    for i = 1 : length(interval(1,:))
        if i == 1
            eje = interval(1,i) : 1 : interval(2,i);
            xim = xo(i) * ones(1,length(eje));
            yim = yo(i) * ones(1,length(eje));
            zim = zo(i) * ones(1,length(eje));
        end
        m = interval(1,i) : 1 : interval(2,i);
        xt=xo(i) * ones(1,length(m));
        yt=yo(i) * ones(1,length(m));
        zt=zo(i) * ones(1,length(m));
        xim = [xim, xt];
        yim = [yim, yt];
        zim = [zim, zt];
        eje = [eje, m];
    end
    ejer = 1 : 1 : length(x);
    
    figure;
    set(gcf,'position',[40 80 1280 720])
    subplot(3,1,1); 
    plot(ejer, x, eje, xim, '.r');   
    legend('X axis', sprintf('Quasi-static periods (%d)', detected));
    xlabel('Sample');
    ylabel('Acceleration (raw)')
    
    subplot(3,1,2); 
    plot(ejer, y, eje, yim, '.r');   
    legend('Y axis', sprintf('Quasi-static periods (%d)', detected));
    xlabel('Sample');
    ylabel('Acceleration (raw)')
    
    subplot(3,1,3); 
    plot(ejer, z, eje, zim, '.r');   
    legend('Z axis', sprintf('Quasi-static periods (%d)', detected));
    xlabel('Sample');
    ylabel('Acceleration (raw)')
    title('DETECTION OF QUASI-STATIC PERIODS')
    
    set(gcf,'NextPlot','add');
    axes;
    h = title('DETECTION OF QUASI-STATIC PERIODS');
    set(gca,'Visible','off');
    set(h,'Visible','on'); 
    
    saveas(gcf, 'figures/calibration/trunk_accRaw_staticDetec.fig')
end
else
    xo = 0;
    yo = 0;
    zo = 0;
end
end
% END OF FUNCTION SELECTSTATICPOSITIONS3D

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [a_parallel, a_antiparallel] = get_acc_parallel_values(acc)

% FUNCTION GET_ACC_PARALLEL_VALUES Extracts the two raw acceleration values
% gathered when the accelerometer axis is placed parallel and antiparallel 
% to the gravity vector. This function requires interaction with the user 
% as it is shown the raw acceleration signal and he has to select the 
% initial and final points of both parallel and anti-parallel static 
% positions. 
% 
% Input parameters:
% |_ 'acc': Raw acceleration signal gathered from the calibration
%           maneuvers.
%
% Output parameters:
% |_ 'a_parallel':      Mode of the gathered acceleration when the 
%                       accelerometer axis is set parallel to the gravity
%                       vector.
% |_ 'a_antiparallel':  Mode of the gathered acceleration when the 
%                       accelerometer axis is set anti-parallel to the 
%                       gravity vector.
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

static_flag = 0;
while static_flag == 0
    
    % 1) Find the indexes of both static periods.
    fig_title = (['SELECT THE INITIAL AND FINAL POINTS OF BOTH THE ',...
        'PARALLEL AND ANTIPARALLEL POSITIONS']);
    indexes = getDCindexes(acc, fig_title);
    
    % 2) Check there are 4 indexes (2 starting points and 2 ending points).
    if length(indexes) == 4
        static_flag = 1;
    else
        msg = msgbox(sprintf(['You have to select 4 points. You ',...
            'selected %d points'],length(indexes)));
        uiwait(msg);
    end
    
end

% 3) Extract the mode of the acceleration value in both parallel and
%    anti-parallel positions. 
bias_chunk1 = mode(acc(indexes(1) : indexes(2)));
bias_chunk2 = mode(acc(indexes(3) : indexes(4)));

% 4) We need to identify the value corresponding to the parallel and
%    antiparallel positions. The value of the parallel position is larger
%    thant the value of the antiparallel position.
if bias_chunk1 > bias_chunk2
    a_parallel = bias_chunk1;
    a_antiparallel = bias_chunk2;
elseif bias_chunk1 < bias_chunk2
    a_parallel = bias_chunk2;
    a_antiparallel = bias_chunk1;
end

end
% END OF GET_ACC_PARALLEL_VALUES FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function euler = quatTOeuler(quatInput)

% FUNCTION FUNCTION QUAT2EULER Computes the Euler angles from a given 
% quaternion.
%
% Input parameters:
% |_ 'quatInput': Nx4 matrix containing a quaternion per row. Each row will 
%                 be converted to a 1x3 matrix containing Euler's 
%                 corresponding angles.
% Output parameters:
% |_'euler': Nx3 matrix containing corresponding Euler's angles for each
%            input quaternion. Angles ar as follows:
%               - euler(:,1): rotation about the X-axis (roll).
%               - euler(:,2): rotation about the Y-axis (pitch).
%               - euler(:,3): rotation about the Z-axis (yaw).
%
% -------------------------------------------------------------------------
% Authors:  Alberto Olivares & Gonzalo Ruiz & Kai B?tzel.
% Entity:  Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------
%
% 1) Check input parameters:
if size(quatInput, 2) ~= 4
    error('Input parameter ''quatInput'' must be a Nx4 matrix.');
end

% 2) Initialize variables.
len = size(quatInput, 1);
euler = zeros(len, 3);

% 3) Apply the transformation equations.
for i=1:1:len
    q0 = quatInput(i, 1);
    q1 = quatInput(i, 2);
    q2 = quatInput(i, 3);
    q3 = quatInput(i, 4);
    euler(i, 1) = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1*q1 + q2*q2));
    euler(i, 2) = asin(2*(q0*q2 - q3*q1));
    euler(i, 3) = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2*q2 + q3*q3));
end

end
% END OF FUNCTION QUATTOEULER

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function quat = eulerToQuat(roll, pitch, yaw)

% FUNCTION EULERTOQUAT computes the orientation quaternion given the three
% euler rotation angles.
%
% Input parameters:
% |_ 'roll':  Rotation angle around the X axis.
% |_ 'pitch': Rotation angle around the Y axis.
% |_ 'yaw':   Rotation angle around the Z axis.
% 
% Output parameters:
% |_ 'quat': Orientation quaternion.
%
% -------------------------------------------------------------------------
% Authors:  Alberto Olivares & Gonzalo Ruiz & Kai B?tzel.
% Entity:  Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Apply transformation equations.
c1 = cos(yaw/2);
s1 = sin(yaw/2);
c2 = cos(pitch/2);
s2 = sin(pitch/2);
c3 = cos(roll/2);
s3 = sin(roll/2);
c1c2 = c1*c2;
s1s2 = s1*s2;
w = c1c2*c3 - s1s2*s3;
x = c1c2*s3 + s1s2*c3;
y = s1*c2*c3 + c1*s2*s3;
z = c1*s2*c3 - s1*c2*s3;

% 2) Built quaternion.
quat = [w x y z];

end
% END OF EULERTOQUAT FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function quat = quat_mul(q1, q2)

% FUNCTION QUAT_MUL Carries out the multiplication between two given
% quaternions.
% 
% Input parameters:
% |_ q1: a 1x4 matrix representing the first quaternion.
% |_ q2: a 1x4 matrix representing the second quaternion.
%
% Output parameters:
%  |_ quat: quaternion resulting from q1 x q2.
%
% -------------------------------------------------------------------------
% Authors:  Alberto Olivares & Gonzalo Ruiz & Kai B?tzel.
% Entity:  Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Check input parameters.
if size(q1, 1)~=1 || size(q1, 2)~=4
    error('Input parameter ''q1'' must be a 1x4 matrix.');
end
if size(q2, 1)~=1 || size(q2, 2)~=4
    error('Input parameter ''q2'' must be a 1x4 matrix.');
end

% 2) Initialize variables.
quat = zeros(1, 4);

% 3) Computation of Hamiltonian product.
quat(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
quat(2) = q1(2)*q2(1) + q1(1)*q2(2) - q1(4)*q2(3) + q1(3)*q2(4);
quat(3) = q1(3)*q2(1) + q1(4)*q2(2) + q1(1)*q2(3) - q1(2)*q2(4);
quat(4) = q1(4)*q2(1) - q1(3)*q2(2) + q1(2)*q2(3) + q1(1)*q2(4);

end
% END OF QUAT_MUL FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [] = plotData(showPlot,x,y,dir_path,physMag,physUnits,segment,...
    plotNames)

% FUNCTION PLOTDATA Plots X vs. Y. It sets the labels and legends specified
% by the user and saves them in three formats (.fig, .png and .eps) in the
% directory specified by the user.
%
% Input parameters:
% |_ 'showPlot':  String containing a flag to indicate whether the figures
%                 are displayed or not. Its possible values are 'on' and 
%                 'off'.
% |_ 'x':         Data vector plotted in the X axis.
% |_ 'y':         Data vector plotted in the Y axis.
% |_ 'dir_path':  Directory where the figures are stored.
% |_ 'physMag':   Name of the physical magnitude plotted in the Y axis.
% |_ 'physUnits': Units of the physical magnitude plotted in the Y axis.
% |_ 'segment':   Segment of the body where the data have been monitored.
% |_ 'plotNames': Legend of the plots displayed in the figure.
% 
% Output parameters: 
% |_ This function returns no parameters. 
%
% -------------------------------------------------------------------------
% Author: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Check input parameters:
if length(x)~=length(y)
    error('X axis and Y axis vectors must be the same length.');
end
if strcmp(showPlot,'on')~=1 && strcmp(showPlot,'off')~=1
    error('The ''showPlot'' flag must be either ''on'' or ''off''.');
end

% 2) Generate figure with provided signals.
h=figure('Visible',showPlot);
hold on;
size_y=size(y);
colors={'b','r','g','black','cyan'};
for i=1:size_y(2);
    g(i)=plot(x,y(:,i),colors{i},'linewidth',2,'DisplayName',plotNames(i));
    h_legend =legend('-DynamicLegend');
end
set(h_legend,'FontSize',18);
set(gca,'FontSize',18);
xlabel('Time (s)','fontsize',24);
ylabel_string=strcat(physMag,' ',physUnits);
ylabel(ylabel_string,'fontsize',24);

% 3) Build figure names and save them in different formats.
if ispc
    completePath=strcat(dir_path,'\',physMag,{' '},segment);
else
    completePath=strcat(dir_path,'/',physMag,{' '},segment);
end
saveas(h,completePath{1},'fig');
saveas(h,completePath{1},'png');
if ispc
    completePath=strcat(dir_path,'\',physMag,{' '},segment,'.eps');
else
    completePath=strcat(dir_path,'/',physMag,{' '},segment,'.eps');
end
saveas(h,completePath{1},'epsc2');

end
% END OF PLOTDATA FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **                ORIENTATION ESTIMATION FUNCTIONS                     **
% *************************************************************************

function [roll, pitch, yaw, q_output] = quat9dofEKF(ax, ay, az, gx, gy,...
    gz, hx, hy, hz, gyroVarX, gyroVarY, gyroVarZ, alpha, mu_gain, f,...
    state_observed)

% FUNCTION QUAT9DOFKF computes a quaternion based 9DOF Kalman filter to 
% estimate the orientation quaternion and the associated Euler angles 
% (roll, pitch and yaw). 
%
% Input parameters:
% |_ 'ax':          Acceleration measured along X axis (N x 1 vector).
% |_ 'ay':          Acceleration measured along Y axis (N x 1 vector).
% |_ 'az':          Acceleration measured along Z axis (N x 1 vector).
% |_ 'gx':          Angular rate measured along X axis (N x 1 vector).
% |_ 'gy':          Angular rate measured along Y axis (N x 1 vector).
% |_ 'gz':          Angular rate measured along Z axis (N x 1 vector).
% |_ 'hx':          Magnetic field measured along X axis (N x 1 vector).
% |_ 'hy':          Magnetic field measured along Y axis (N x 1 vector).
% |_ 'hz':          Magnetic field measured along Z axis (N x 1 vector).
% |_ 'gyroVarX':    Variance of X-axis gyroscope noise (scalar).
% |_ 'gyroVarY':    Variance of Y-axis gyroscope noise (scalar).
% |_ 'gyroVarZ':    Variance of Z-axis gyroscope noise (scalar).
% |_ 'alpha':       Variance of measurement noise (scalar).
% |_ 'mu_gain':     Gain of the gradient descent step length (scalar).
% |_ 'f':           Sampling frequency (scalar).
% |_ 'state_observed': Initial value of the observed orientation 
%                   quaternion.
% Output parameters:
% |_ 'roll':     Roll angle (N x 1 vector).
% |_ 'pitch':    Pitch angle (N x 1 vector).
% |_ 'yaw':      Yaw angle (N x 1 vector).
% |_ 'q_output': Orientation quaternion (N x 4 matrix).
% 
% -------------------------------------------------------------------------
% Authors:  Alberto Olivares & Kai B?tzel (based on the Visual Studio 
%           version from the Microelectronics Laboratory of University 
%           of Bergamo).
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Definition and initialization of variables.

dt = 1/f;       % Sampling period. 

gx_old = 0;     % Initial value of angular rate around X axis.
gy_old = 0;     % Initial value of angular rate around Y axis.
gz_old = 0;     % Initial value of angular rate around Z axis.

H = eye(4);     % Measurement matrix (observation model which maps the true
                % state space into the observed space).
F = zeros(4,4); % State transition matrix.
R = eye(4);     % Measurement noise covariance matrix.
Q = zeros(4,4); % Process noise covariance matrix.

state_filtered = zeros(4,1);  % Initial value of the filtered state vector.

% 2) Transformation of gyro variances into radians.
sigmaRoll = (gyroVarX/180*pi)^2;  
sigmaPitch = (gyroVarY/180*pi)^2;
sigmaYaw = (gyroVarZ/180*pi)^2;

% 3) Population of process noise covariance matrix.
Q(1, 1) = sigmaRoll + sigmaPitch + sigmaYaw;
Q(1, 2) = -sigmaRoll + sigmaPitch - sigmaYaw;
Q(1, 3) = -sigmaRoll - sigmaPitch + sigmaYaw;
Q(1, 4) = sigmaRoll - sigmaPitch - sigmaYaw;
Q(2, 1) = -sigmaRoll + sigmaPitch - sigmaYaw;
Q(2, 2) = sigmaRoll + sigmaPitch + sigmaYaw;
Q(2, 3) = sigmaRoll - sigmaPitch - sigmaYaw;
Q(2, 4) = -sigmaRoll - sigmaPitch + sigmaYaw;
Q(3, 1) = -sigmaRoll - sigmaPitch + sigmaYaw;
Q(3, 2) = sigmaRoll - sigmaPitch - sigmaYaw;
Q(3, 3) = sigmaRoll + sigmaPitch + sigmaYaw;
Q(3, 4) = -sigmaRoll + sigmaPitch - sigmaYaw;
Q(4, 1) = sigmaRoll - sigmaPitch - sigmaYaw;
Q(4, 2) = -sigmaRoll - sigmaPitch + sigmaYaw;
Q(4, 3) = -sigmaRoll + sigmaPitch - sigmaYaw;
Q(4, 4) = sigmaRoll + sigmaPitch + sigmaYaw;

% 4) Population of measurement noise covariance matrix.
R = alpha*R;

% 5) Initialize a posteriori estimate covariance matrix.
P_update = 0.1*eye(4);

% 6) Begin filter loop.
for i = 1: length(ax)
    
    % Read measurements.
    ax_meas = ax(i);
    ay_meas = ay(i);
    az_meas = az(i);
    gx_meas = gx(i);
    gy_meas = gy(i);
    gz_meas = gz(i);
    hx_meas = hx(i);
    hy_meas = hy(i);
    hz_meas = hz(i);
   
    % Normalize measurements.
    norm = sqrt(ax_meas.^2+ay_meas.^2+az_meas.^2);
    ax_meas = ax_meas/norm;
    ay_meas = ay_meas/norm;
    az_meas = az_meas/norm;

    norm = sqrt(hx_meas.^2+hy_meas.^2+hz_meas.^2);
    hx_meas = hx_meas/norm;
    hy_meas = hy_meas/norm;
    hz_meas = hz_meas/norm;
    
    gx_meas = gx_meas/180*pi;
    gy_meas = gy_meas/180*pi;
    gz_meas = gz_meas/180*pi;
    
    % Average the current angular rate measurements with the previous
    % instant measurements.
    gx_meas = (gx_meas + gx_old) / 2;
    gy_meas = (gy_meas + gy_old) / 2;
    gz_meas = (gz_meas + gz_old) / 2;
    
    % Compute and normalize the quaternion derivative: 
    %  |_ 'dq = 1/2*(q_obs x q_g)'.
    % dq = 0.5 * quatmultiply(state_observed',[0 gx_meas gy_meas gz_meas]);
    dq = 0.5 * quat_mul(state_observed',[0 gx_meas gy_meas gz_meas]);
    norm = sqrt(dq(1)^2 + dq(2)^2 + dq(3)^2 + dq(4)^2);

    % Compute the gradient descent step size.
    mu = mu_gain * norm * dt;
    
    % Compute state observation (Accelerometer + Magnetometer).
    state_observed = gradientDescent(ax_meas, ay_meas, az_meas, hx_meas,...
        hy_meas, hz_meas, mu, state_observed); 

    % Populate state transition matrix.
    F(1, 1) = 1;
    F(1, 2) = -dt / 2 * gx_meas;
    F(1, 3) = -dt / 2 * gy_meas;
    F(1, 4) = -dt / 2 * gz_meas;
    F(2, 1) = dt / 2 * gx_meas;
    F(2, 2) = 1;
    F(2, 3) = dt / 2 * gz_meas;
    F(2, 4) = -dt / 2 * gy_meas;
    F(3, 1) = dt / 2 * gy_meas;
    F(3, 2) = -dt / 2 * gz_meas;
    F(3, 3) = 1;
    F(3, 4) = dt / 2 * gx_meas;
    F(4, 1) = dt / 2 * gz_meas;
    F(4, 2) = dt / 2 * gy_meas;
    F(4, 3) = -dt / 2 * gx_meas;
    F(4, 4) = 1;

    % Compute state prediction using the filtered state vector from
    % previous step.
    state_predicted = F * state_filtered;

    % Compute the a priori estimate covariance matrix.
    P_predicted = F * P_update * F' + Q;

    % Compute Kalman gain.
    temp = H * P_predicted * H' + R;
    K = P_predicted * H' * temp;

    % Update the state vector of the system.
    temp = state_observed - H * state_predicted;
    state_filtered = state_predicted + K * temp;

    % Compute the a posteriori estimate covariance matrix.
    temp = eye(4) - K * H;
    P_update = temp * P_predicted;

    % Normalize state vector.
    norm = sqrt(state_filtered(1)^2 + state_filtered(2)^2 + ...
        state_filtered(3)^2 + state_filtered(4)^2);
    state_filtered = state_filtered / norm; 
    q_output(i,:) = state_filtered';

    % Update the value of the angular rate measurements from previous step.
    gx_old = gx_meas;
    gy_old = gy_meas;
    gz_old = gz_meas;

end

% 7) Transform computed quaterniones into Euler angles.
roll = zeros(1,length(q_output));
pitch = zeros(1,length(q_output));
yaw = zeros(1,length(q_output));

for i = 1: length(q_output)
    eulerAngles = quatTOeuler(q_output(i,:));
    roll(i) = eulerAngles(1);
    pitch(i) = eulerAngles(2);
    yaw(i) = eulerAngles(3);
end

end
% END OF QUAT9DOFKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function result = gradientDescent(ax, ay, az, hx, hy, hz, mu, ...
    state_observed)

% FUNCTION GRADIENTDESCENT computes the observation orientation quaternion 
% using the measured acceleration and magnetic field in the sensor frame. 
% The optimal quaternion is found minimizing an error function which by 
% means of the gradient descent method.
% 
% Iput parameters:
% |_ 'ax': Acceleration measured along X axis (vector).
% |_ 'ay': Acceleration measured along Y axis (vector).
% |_ 'az': Acceleration measured along Z axis (vector).
% |_ 'hx': Magnetic field measured along X axis (vector).
% |_ 'hy': Magnetic field measured along Y axis (vector).
% |_ 'hz': Magnetic field measured along Z axis (vector).
% |_ 'mu': Step size of the gradient descent (scalar).
% |_ 'state_observed': Initial orientation quaternion (4 x 1 vector).
%
% Output parameters:
% |_ 'result': computed orientation quaternion (4 x 1 vector). 
%
% -------------------------------------------------------------------------
% Authors:  Alberto Olivares & Kai B?tzel (based on the Visual Studio 
%           version from the Microelectronics Laboratory of University 
%           of Bergamo).
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 28/11/2013.
% -------------------------------------------------------------------------

% 1) Initialize the value of the objective function.
f_obb = zeros(6,1);

% 2) Initialize the Jacobian matrix.
Jacobian = zeros(6,4);

% 3) Get data from the initial orientation quaternion.
q1 = state_observed(1);
q2 = state_observed(2);
q3 = state_observed(3);
q4 = state_observed(4);

% 4) Build initial value of the orientation quaternion.
quaternion = [q1 q2 q3 q4];
quaternion_conjugate = [q1 -q2 -q3 -q4];
% temp = quatmultiply(quaternion,[0 hx hy hz]);
temp = quat_mul(quaternion,[0 hx hy hz]);
% h = quatmultiply(temp,quaternion_conjugate);
h = quat_mul(temp,quaternion_conjugate);

% 5) Compute the vertical and horizontal components of the magnetic field.
bx = sqrt(h(2)^2 + h(3)^2);
by = 0;
bz = h(4);

norm = sqrt(bx^2 +bz^2);
bx = bx/norm;
bz = bz/norm;

% 6) Compute the objective function.
f_obb(1) = 2 * (q2 * q4 - q1 * q3) - ax;
f_obb(2) = 2 * (q1 * q2 + q3 * q4) - ay;
f_obb(3) = 2 * (0.5 - q2 * q2 - q3 * q3) - az;
f_obb(4) = 2 * bx * (0.5 - q3 * q3 - q4 * q4) + 2 * by * (q1 * q4 + q2 *...
    q3) + 2 * bz * (q2 * q4 - q1 * q3) - hx;
f_obb(5) = 2 * bx * (q2 * q3 - q1 * q4) + 2 * by * (0.5 - q2 * q2 - q4 *...
    q4) + 2 * bz * (q1 * q2 + q3 * q4) - hy;
f_obb(6) = 2 * bx * (q1 * q3 + q2 * q4) + 2 * by * (q3 * q4 - q1 * q2) +...
    2 * bz * (0.5 - q2 * q2 - q3 * q3) - hz;

% 7) Compute the Jacobian matrix.
Jacobian(1, 1) = -2 * q3;
Jacobian(1, 2) = 2 * q4;
Jacobian(1, 3) = -2 * q1;
Jacobian(1, 4) = 2 * q2;
Jacobian(2, 1) = 2 * q2;
Jacobian(2, 2) = 2 * q1;
Jacobian(2, 3) = 2 * q4;
Jacobian(2, 4) = 2 * q3;
Jacobian(3, 1) = 0;
Jacobian(3, 2) = -4 * q2;
Jacobian(3, 3) = -4 * q3;
Jacobian(3, 4) = 0;
Jacobian(4, 1) = 2 * by * q4 - 2 * bz * q3;
Jacobian(4, 2) = 2 * by * q3 + 2 * bz * q4;
Jacobian(4, 3) = -4 * bx * q3 + 2 * by * q2 - 2 * bz * q1;
Jacobian(4, 4) = -4 * bx * q4 + 2 * by * q1 + 2 * bz * q2;
Jacobian(5, 1) = -2 * bx * q4 + 2 * bz * q2;
Jacobian(5, 2) = 2 * bx * q3 - 4 * by * q2 + 2 * bz * q1;
Jacobian(5, 3) = 2 * bx * q2 + 2 * bz * q4;
Jacobian(5, 4) = -2 * bx * q1 - 4 * by * q4 + 2 * bz * q3;
Jacobian(6, 1) = 2 * bx * q3 - 2 * by * q2;
Jacobian(6, 2) = 2 * bx * q4 - 2 * by * q1 - 4 * bz * q2;
Jacobian(6, 3) = 2 * bx * q1 + 2 * by * q4 - 4 * bz * q3;
Jacobian(6, 4) = 2 * bx * q2 + 2 * by * q3;

% 8) Compute gradient of the function.
Df = Jacobian'*f_obb;
norm = sqrt(Df(1)^2 + Df(2)^2 + Df(3)^2 + Df(4)^2);

% 9) Compute norm of the gradient.
Df = Df / norm;

% 10) Update the orientation quaternion using the gradient descent.
result = quaternion' - mu * Df;

% 11) Renormalize the orientation quaternion.
norm = sqrt(result(1)^2 + result(2)^2 + result(3)^2 + result(4)^2);
result = result / norm;

end
% END OF FUNCTION GRADIENTDESCENT

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimize_KF(gyro, acc, frec, obsVar,...
    accVar, gyroVar, Ini, ref_angle, p0, rmse_off)

% FUNCTION OPTIMIZEKF uses the ANMS algorithm to find the optimal 
% parameters of the Regular Kalman Filter which minimize the error between
% the actual orientation angle and the one estimated with KF. 
% 
% - INPUT PARAMETERS:
%   |_ 'observation': Observation of the Kalman Filter.
%   |_ 'ang_vel': Angular velocity (measured with the gyroscope).
%   |_ 'frequency': Sampling frequency.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'rmse_off': Offset in the RMSE computation.
%
% - OUPUT PARAMETERS: 
%   |_ 'xmin': Value of the parameters which minimize the error function.
%   |_ 'fmin': Minimum value of the error function (minimum RMSE).+
%   |_ 'ct': Number of algorithm iterations to find the minimum.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHOR:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares. All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
%  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------
%
% 1) Set variables.
% -------------------------------------------------------------------------
global gyrosc;        gyrosc = gyro;
global obs;           obs = acc;
global true_angle;    true_angle = ref_angle;
global freq;          freq = frec;
global rmse_offset;   rmse_offset = rmse_off;
global obsVariance;   obsVariance = obsVar;
global accVariance;   accVariance = accVar;
global gyroVariance;  gyroVariance = gyroVar;
global Initial;       Initial = Ini;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp('Optimizing parameters of Kalman Filter (it may take a while) ...');

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin, fmin, ct] = ANMS(@eofKalman, p0, tol, max_feval);

end
% END OF OPTIMIZE_KF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofKalman(p)

% FUNCTION EOFKALMAN Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the Kalman Filter. 
%
% - INPUT PARAMETERS:
%   |_ 'p': Vector of initial value of the parameters to be optimized.
%
% - OUTPUT PARAMETERS:
%   |_ 'F': Value of the error function.
% 
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHOR:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares. All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
%  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------

% 1) Set variables.
% -------------------------------------------------------------------------
global gyrosc;
global obs;
global freq;
global true_angle;
global rmse_offset;
global obsVariance;
global accVariance;
global gyroVariance;
global Initial;

% 2) Extract the first parameter to be minimized (measurement noise 
% variance gain). 
alpha = p(1); 
beta = p(2);

% 3) Estimate the orientation angle using the Kalman Filter.
angle_KF = fusion_KF(gyrosc, obs, freq, obsVariance, ...
           accVariance, gyroVariance, alpha, beta, Initial);

% 4) Compute the error function.
F = sqrt(mean((true_angle(rmse_offset : end) -  ...
    angle_KF(rmse_offset : end)') .^ 2));
end
% END OF EOFKALMAN FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [x1,x2,v_output,K_output] = fusion_KF(gyro, acc, frec, obsVar,...
    accVar, gyroVar, alfa, beta, Ini)

% FUNCTION FUSION_KF Applies a Kalman filter to fuse the accelerometer and
% gyroscope data to obtain an accurate estimate of the position angle.
%
% Input parameters:
% |_ 'gyro':    Vector containing the gyroscope signal for a determined 
%               axis. 
% |_ 'acc':     Vector containing the angle values obtained via the linear
%               acceleration relations.
% |_ 'frec':    Sampling frecuency. Must be real positive.
% |_ 'obsVar':  Variance of the observation/measurement process in the 
%               Kalman's filtering model. Must be real positive.
% |_ 'accVar':  Variance of the acceleration sensor in static conditions. 
%               Must be real positive.
% |_ 'gyroVar': Variance of the gyroscope sensor in static conditions. Must
%               be real positive.
% |_ 'alfa':    Scalar weighing the variance of the observation. This
%               parameter can be adjusted to tune the filter.
% |_ 'beta':    Scalar weighing the variance of the angular rate. This
%               parameter can be adjusted to tune the filter.
%
% Output parameters:
% |_ 'x1':       First element of the estimated state (angle value 
%                corrected by the Kalman filter).
% |_ 'x2':       Second element of the estimated state (gyroscope bias).
% |_ 'v_output': Estimation error of each iteration.
% |_ 'K_output': Kalman gain of each iteration. 

% -------------------------------------------------------------------------
% Authors: Alberto Olivares, Gonzalo Ruiz & Kai B?otzel.
% Entity: Universidad de Granada & Ludwig-Maixmilians Universit?t M?nchen.
% Last modification: 29/05/2012 .
% -------------------------------------------------------------------------

% 1) Check input parameters:
len = length(gyro);
if length(acc)~=len
    error(['Input parameters ''gyro'' and ''acc'' must be the same ' ...
        'length.']);
end
if frec <= 0
    error('Input parameter ''frec'' must be real positive.');
end
if obsVar <= 0
    error('Input parameter ''obsVar'' must be real positive.');
end
if accVar <= 0
    error('Input parameter ''accVar'' must be real positive.');
end
if gyroVar <= 0
    error('Input parameter ''gyroVar'' must be real positive.');
end


% 1) Definition of the measurement noise covariance and the sampling
%    period.
Rk = alfa*obsVar;
dt = 1/frec;

% 2) Initialization of covariance matrix.
P = [1 0; 0 1];

% 3) Assuming process noise is white for each component, covariance matrix 
%    of process noise must be a diagonal matrix, with noise variance for 
%    each component in each position of the diagonal.
Q = [beta*gyroVar 0; 0 beta*accVar];

% 4) Definition of state transition matrix.
A = [1 -dt; 0 1];

% 5) Definition of the measurement matrix (which relates the measurements 
%    to the state vector).
C = [1 0];

% 6) Initializing state vector.
X = [Ini; gyro(1)];
x1 = zeros(len, 1);
x2 = zeros(len, 1);

% 7) Begin filter loop.
for i=1:1:len-1
        
        % Prediction phase
        X(1) = X(1) + (gyro(i)-X(2))*dt;
        P = A*P*A' + Q;

        % Update phase
        v = acc(i) - C*X;
        v_output(i) = v;
        Sk = C*P*C' + Rk;
        K = (P*C')/Sk;
        K_output(i,:) = K';
        X = X + K*v;
        P = P - K*Sk*K';
        x1(i) = X(1);
        x2(i) = X(2);   
end

x1(end) = x1(end-1);

end
% END OF FUSION_KF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [x_min, f_min, ct] = optimize_EKF( ...
                        gyro_thigh_y, gyro_shank_y, ...
                        acc_thigh_x, acc_thigh_z, ...
                        acc_shank_x, acc_shank_z, ...
                        fs, l1, l2, ref_angles, p0, ...
                        rmse_off)

% FUNCTION OPTIMIZE_EKF uses an adaptive Nelder-Mead 
% simplex ANMS algorithm to find the optimal parameters
% of the Extended Kalman filter. 
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
% |_ 'p':             Row vector consisting of the 
%                     filter parameters sigma_t1, 
%                     sigma_t2, sigma_b, sigma_f_1, 
%                     sigma_f_2, sigma_s_1, sigma_s_2.
% |_ 'ref_angle':     Orientation angle reference.
% |_ 'p0':            Initial value of the parameters 
%                     to be optimized.
% |_ 'rmse_off':      Offset in the RMSE computation.
%
% % Output:
% |_ 'xmin':          Value of the parameters that 
%                     minimize the error function.
% |_ 'fmin':          Minimum value of the error 
%                     function.
% |_ 'ct':            Number of algorithm iterations
%                     to find the minimum.

% 1) Set variables.
% -----------------------------------------------------
global gyro_thigh_y_g; gyro_thigh_y_g = gyro_thigh_y;
global gyro_shank_y_g; gyro_shank_y_g = gyro_shank_y;
global acc_thigh_x_g; acc_thigh_x_g = acc_thigh_x;
global acc_thigh_z_g; acc_thigh_z_g = acc_thigh_z;
global acc_shank_x_g; acc_shank_x_g = acc_shank_x;
global acc_shank_z_g; acc_shank_z_g = acc_shank_z;
global fs_g; fs_g = fs;
global l1_g; l1_g = l1;
global l2_g; l2_g = l2;
global true_angles;    true_angles = ref_angles;
global rmse_offset;   rmse_offset = rmse_off;

% 2) Call the minimisation routine.
% -----------------------------------------------------
disp('Optimising parameters of extended Kalman Filter...');

% Set tolerance limit between the minimum values found 
% in subsequent iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error
% function.
max_feval = 5000;

% Call ANMS algorithm which minimises the error function.
[x_min, f_min, ct] = ANMS(@eofEKF, p0, tol, max_feval);

end
% END OF OPTIMIZE_EKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofEKF(p)

% FUNCTION EOFKALMAN is the error function to be 
% minimised: That is, the sum of the RMSE between the 
% actual angle and the estiamted orientation angle 
% computed with the extended Kalman filter. 
%
% Input arguments:
% |_ 'p':   Vector of initial value of the parameters 
%           to be optimized.
%
% Output:
% |_ 'F':   Value of the error function.

% 1) Set variables.
% -----------------------------------------------------
global gyro_thigh_y_g; 
global gyro_shank_y_g;
global acc_thigh_x_g;
global acc_thigh_z_g;
global acc_shank_x_g;
global acc_shank_z_g;
global fs_g;
global l1_g;
global l2_g;
global true_angles;
global rmse_offset;

% 3) Estimate the orientation angle using the extended 
%    Kalman Filter.
[theta1, theta2, ~, ~, ~] = fusion_EKF(...
                   gyro_thigh_y_g, gyro_shank_y_g, ...
                   acc_thigh_x_g, acc_thigh_z_g, ...
                   acc_shank_x_g, acc_shank_z_g, ...
                   fs_g, l1_g, l2_g, p);
               
thigh_angle_EKF = theta1;
shank_angle_EKF = theta1 + theta2;

% 4) Compute the error function.
F1 = sqrt(mean((true_angles(1, rmse_offset : end) -  ...
     thigh_angle_EKF(rmse_offset : end)) .^ 2));
F2 = sqrt(mean((true_angles(2, rmse_offset : end) -  ...
     shank_angle_EKF(rmse_offset : end)) .^ 2));

F = F1 + F2;
end
% END OF EOEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


function corrected_angle = correct_quad_shifts(angle,units)

% FUNCTION CORRECT_QUAD_SHIFTS corrects shifts of ?360? (or ?2pi rad) in 
% Euler angle signals. 
%
% Input parameters:
% |_ 'angle': Signal containing the Euler angles (vector).
% |_ 'units': Units of the angle (degrees or radians).
%
% Ouput parameters:
% |_ 'corrected_angle': Corrected angle signal without quadran shifts
%                       (vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 07/11/2013.
% -------------------------------------------------------------------------

% 1) Input parameters.
if ~(strcmpi(units,'deg') || strcmpi(units,'rad'))
  disp('ERROR: Unites must be either degrees ("deg") or radians ("rad").');
end

% 2) Define the correction angle shift.
corrected_angle = angle;
if strcmpi(units,'deg')
    angle_shift = 360;
elseif strcmpi(units,'rad')
    angle_shift = 2*pi;
end

% 3) Check if the angles need to be corrected.
for i = 1:length(angle)
    if angle(i) < - 290
        corrected_angle(i) = angle(i) + angle_shift; 
    end
end

end
% END OF CORRECT_QUAD_SHIFTS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function corrected_angle = correct_yaw_quad_shifts(uncorrected_angle,units)

% CORRECT_QUAD_SHIFTS corrects shifts of ?360? (or ?2pi rad) in Euler
% angle signals. 
%
% Input parameters:
% |_ 'angle': Signal containing the Euler angles (vector).
% |_ 'units': Units of the angle (degrees or radians).
%
% Ouput parameters:
% |_ 'corrected_angle': Corrected angle signal without quadran shifts
%                       (vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 22/01/2014.
% -------------------------------------------------------------------------

% Check angle units.
if ~(strcmpi(units,'deg') || strcmpi(units,'rad'))
  disp('ERROR: Unites must be either degrees ("deg") or radians ("rad").');
end

% Define the minimum angle shift between consecutive angles to define a
% quadrant shift, and the correction angle. This is done for both degrees
% and radians.
% if strcmpi(units,'deg')
%     limit = 356;
%     correction_angle = 360;
%     final_correction = 90;
% elseif strcmpi(units,'rad')
%     limit = 356/180*pi;
%     correction_angle = 2*pi;
%     final_correction = pi/2;
% end

% % Map all angles to [0,360)
% for l=1:length(uncorrected_angle)
%     if uncorrected_angle(l) < 0
%         uncorrected_angle(l) = uncorrected_angle(l) + correction_angle;
%     end
% end
% 
% corrected_angle = uncorrected_angle;
% indexes = [];
% 
% % The signal is swept and we look for quadrant shifts. The indexes in which
% % the quadrant shift starts and finishes are stored. 
% for i = 1: length(uncorrected_angle)-1
%     if abs(uncorrected_angle(i+1) - uncorrected_angle(i)) >= limit
%         indexes = [indexes i];
%     end
% end
% 
% % Correct the quadrant shifts. 
% for i = 1:2: length(indexes)-1
%     if uncorrected_angle(indexes(i)+1) > limit 
%         corrected_angle(indexes(i)+1:indexes(i+1)) = uncorrected_angle(indexes(i) + ...
%             1:indexes(i+1)) - correction_angle;
%     elseif uncorrected_angle(indexes(i)+1) < - limit
%         corrected_angle(indexes(i)+1:indexes(i+1)) = uncorrected_angle(indexes(i) + ...
%             1:indexes(i+1)) + correction_angle;
%     end
% end
% 
% % if the number of indexes is odd, that means that the last quadrant shift
% % does not have an ending point (and therefore the ending point will be the
% % last sample of the angle signal).
% if mod(length(indexes),2) ~= 0
%     if uncorrected_angle(indexes(end)+1) > limit
%         corrected_angle(indexes(end)+1:end) = uncorrected_angle(indexes(end)+1:end) -...
%             correction_angle;
%     elseif uncorrected_angle(indexes(end)+1) < - limit
%         corrected_angle(indexes(end)+1:end) = uncorrected_angle(indexes(end)+1:end) +...
%             correction_angle;
%     end
% end
% 
% corrected_angle = corrected_angle - final_correction;
limit = 250;
correction = 0;
if min(uncorrected_angle) < 0
    corrected_angle = uncorrected_angle + abs(min(uncorrected_angle));
    correction = 1;
else
    corrected_angle = uncorrected_angle;
end

for i = 2: length(corrected_angle)
    if abs(corrected_angle(i) - corrected_angle(i-1)) > limit
        if corrected_angle(i) > corrected_angle(i-1) 
            corrected_angle(i) = corrected_angle(i) - 360;
        elseif corrected_angle(i) < corrected_angle(i-1) 
            corrected_angle(i) = corrected_angle(i) + 360;
        end
    end
end
if correction == 1
    corrected_angle = corrected_angle - abs(min(uncorrected_angle));
end
end
% END OF CORRECT_QUAD_SHIFTS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function corrected_angle = correct_ekf_quad_shifts(uncorrected_angle,units)

% CORRECT_QUAD_SHIFTS corrects shifts of ?360? (or ?2pi rad) in Euler
% angle signals. 
%
% Input parameters:
% |_ 'angle': Signal containing the Euler angles (vector).
% |_ 'units': Units of the angle (degrees or radians).
%
% Ouput parameters:
% |_ 'corrected_angle': Corrected angle signal without quadran shifts
%                       (vector).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 22/01/2014.
% -------------------------------------------------------------------------

% Check angle units.
if ~(strcmpi(units,'deg') || strcmpi(units,'rad'))
  disp('ERROR: Unites must be either degrees ("deg") or radians ("rad").');
end

% Define the minimum angle shift between consecutive angles to define a
% quadrant shift, and the correction angle. This is done for both degrees
% and radians.
if strcmpi(units,'deg')
    limit = 350;
    correction_angle = 360;
elseif strcmpi(units,'rad')
    limit = 350/180*pi;
    correction_angle = 2*pi;
end

corrected_angle = uncorrected_angle;
indexes = [];

% The signal is swept and we look for quadrant shifts. The indexes in which
% the quadrant shift starts and finishes are stored. 
for i = 1: length(uncorrected_angle)-1
    if abs(uncorrected_angle(i+1) - uncorrected_angle(i)) >= limit
        indexes = [indexes i];
    end
end

% Correct the quadrant shifts. 
for i = 1:2: length(indexes)-1
    if uncorrected_angle(indexes(i)+1) > limit
        corrected_angle(indexes(i)+1:indexes(i+1)) = uncorrected_angle(indexes(i) + ...
            1:indexes(i+1)) - correction_angle;
    elseif uncorrected_angle(indexes(i)+1) < - limit
        corrected_angle(indexes(i)+1:indexes(i+1)) = uncorrected_angle(indexes(i) + ...
            1:indexes(i+1)) + correction_angle;
    end
end

% if the number of indexes is odd, that means that the last quadrant shift
% does not have an ending point (and therefore the ending point will be the
% last sample of the angle signal).
if mod(length(indexes),2) ~= 0
    if uncorrected_angle(indexes(end)+1) > limit
        corrected_angle(indexes(end)+1:end) = uncorrected_angle(indexes(end)+1:end) -...
            correction_angle;
    elseif uncorrected_angle(indexes(end)+1) < - limit
        corrected_angle(indexes(end)+1:end) = uncorrected_angle(indexes(end)+1:end) +...
            correction_angle;
    end
end
end
% END OF CORRECT_QUAD_SHIFTS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function int_rate = integRate(dt,g,ini)

% FUNCTION INTEGRATE Carries out the integration of the values contained in
% a vector with a specific sampling period.
%
% Input parameters:
% |_ 'dt':  Sampling period.
% |_ 'g':   Vector to be integrated.
% |_ 'ini': Initial value of the integrated vector.
%
% Output parameters:
% |_ 'int_rate': Integrated vector.
%
% -------------------------------------------------------------------------
% Author: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 13/12/2013.
% -------------------------------------------------------------------------

% 1) Initialize variables.
int_rate=zeros(1,length(g));
int_rate(1)=ini;  
dt_int=ini;

% 2) Integrate angular rate.
for i=2:length(g)
        dt_int = dt_int + 0.5*dt*(g(i) + g(i-1));
        int_rate(i) = dt_int;
end

end
% END OF INTEGRATE FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% *************************************************************************
% **                            MISCELLANEA                              **
% *************************************************************************

function [index_GW, index_QS] = sync_GW_QS(gw_signal, qs_signal)

% FUNCTION SYNC_GW_QS is used to synchronize GaitWatch data with Qualisys 
% Tracking data. It shows the user the two signals and he has to select the
% first maximum or minimum point. These two indexes can be used to apply a 
% time shift to synchronize the signals. 
%
% Input parameters:
% |_ 'gw_signal': Signal containing the GaitWatch data (vector).
% |_ 'qs_signal': Signal containing the Qualisys data (vector).
%
% Output parameters:
% |_ 'index_GW': Index of first maximum/minimum in GaitWatch signal 
%                (scalar).
% |_ 'index_QS': Index of first maximum/minimum in Qualisys signal 
%                (scalar).
%
% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai B?tzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universit?t M?nchen.
% Last modification: 21/11/2013.
% -------------------------------------------------------------------------

% 1) Select the synchronization point of the GaitWatch signal.
flag1 = 0;
while flag1 == 0
    fig_title = ['Select the first maxima or minima from which you wish',...
    ' to synchronize the signals (GaitWatch)'];
    index_GW = getDCindexes(gw_signal,fig_title);
    if length(index_GW) == 1
        flag1 = 1;
    end
end
close all

% 2) Select the synchronization point of the Qualisys signal.
flag2 = 0;
while flag2 == 0
    fig_title = ['Select the first maxima or minima from which you wish',...
    ' to synchronize the signals (Qualisys)'];
    index_QS = getDCindexes(qs_signal,fig_title);
    if length(index_QS) == 1
        flag2 = 1;
    end
end
close all

end

% *************************************************************************
% **            MOTION INTENSITY DETECTION FUNCTIONS                     **
% *************************************************************************

function SP = spectrum(x, LWIN, SHIFT, NFFT)

% Calculo el espectro de la se?al x para cada frame.

Nframes= floor((length(x)-LWIN)/SHIFT)+1;

SP= zeros(NFFT/2, Nframes);

for frame=1:Nframes
    r=(frame-1)*SHIFT+[1:LWIN];
    s_w= hamming(LWIN).*x(r);
    %s_w= blackman(LWIN).*x(r);    
    magnitude= abs(fft(s_w,NFFT));
    SP(:,frame)= magnitude(1:NFFT/2);
end
end

function [V,T,Threshold,Nframes] = fsd(x,LWIN,SHIFT,NFFT,Threshold)


Nframes= floor((length(x)-LWIN)/SHIFT)+1;
% Computation of spectrum.
SP = spectrum(x,LWIN,SHIFT,NFFT);

V= zeros(Nframes,1);
% Noise averaging initizalization time (frames)
FI= 1;
% Noise estimation
NE=0;
alfa= 0.98;
% Decision function
T= zeros(Nframes,1);

for frame=1:Nframes
       
    if (frame<=FI)
       NE= (1-1/frame)*NE + 1/frame*SP(:,frame);
       V(frame)= 0;
    else
       T(frame)= 10*log10(1/(NFFT/2)*sum(SP(:,frame).^2./NE.^2));
       if (T(frame)>Threshold)
           V(frame)= 1; 
       else
           V(frame)= 0;
           NE= alfa*NE+(1-alfa)*SP(:,frame);
       end
    end
end
end

function [y1,y2] = compEstMark(V, T, x, LWIN, SHIFT)
%--------------------------------------------------------------------------
% Input parameters:
% -V: Figure of merit computed with LTSD or FSD.
% -x: Signal used as input to compute the figure of merit with LTSD or FSD,
% i.e. acceleration or angular rate magnitude, sum of magnitudes or product
% of magnitudes.
% -LWIN: Window length used to compute the figure of merit with LTSD or
% FSD.
% -SHIFT: Overlapping used to computed the figure of merit with LTSD or
% FSD.
% Output parameters:
% -y: (In)activity binary marker.
%--------------------------------------------------------------------------
y1=zeros(1,length(x));
y2=zeros(1,length(x));

N= min(floor((length(x)-LWIN)/SHIFT)+1, length(V));

for i=1:N
    y2((i-1)*SHIFT+(1:LWIN))= T(i);
   if (V(i)==0)
          y1((i-1)*SHIFT+(1:LWIN))= 0; 
   else
          y1((i-1)*SHIFT+(1:LWIN))= 1; 
   end
end
end

function [x1,x2,P_output] = fusionGKF(pGyro, pAcc, pFrec, pObsVar, pAccVar, ...
    pGyroVar,alfa1,alfa2,beta1,beta2,Ini,marker)
%% HELP:
% INPUT PARAMETERS
% -------------------------------------------------------------------------
%   gyro: vector containing the gyroscope signal for a determined axis. 
%   acc: vector containing the angle values obtained via the linear
%       acceleration relations.
%   frec: sampling frecuency. Must be real positive.
%   obsVar: variance of the observation/measurement process in the Kalman's
%       filtering model. Must be real positive.
%   accVar: variance of the acceleration sensor in static conditions. Must
%       be real positive.
%   gyroVar: variance of the gyroscope sensor in static conditions. Must be
%       real positive.
%   angleTransition: value of angle where a discontinuity overcomes due to
%       inverse trigonometric functions. Typical values are \pi or \pi/2.
%       Must be real positive.
%   angleThres: threshold about angleTransition to ignore filtering and
%       avoid Kalman filter to smooth transitions. Must be real positive.
%
%   IMPORTANT NOTE: gyro and acc must be the same length. Otherwise, an
%       error wil be returned.
%   
%
% OUTPUT PARAMETERS
% -------------------------------------------------------------------------
%   x1: angle values corrected by Kalman filter.

% CHECKING INPUT PARAMETERS
% Later
len = length(pGyro);
if length(pAcc)~=len
    error(['Input parameters ''gyro'' and ''acc'' must be the same ' ...
        'length.']);
end
if pFrec <= 0
    error('Input parameter ''frec'' must be real positive.');
end
if pObsVar <= 0
    error('Input parameter ''obsVar'' must be real positive.');
end
if pAccVar <= 0
    error('Input parameter ''accVar'' must be real positive.');
end
if pGyroVar <= 0
    error('Input parameter ''gyroVar'' must be real positive.');
end

% LOCAL COPY OF INPUT PARAMETER
gyro = pGyro;
acc = pAcc;
frec = pFrec;
Rk1 = alfa1*pObsVar;
Rk2 = alfa2*pObsVar;
% Rk = 20*pObsVar;
accVar = pAccVar;
gyroVar = pGyroVar;

% Defining useful variables
dt = 1/frec;

% State vector will be x = [angle; angle_bias].

% Initialization of covariance matrix.
P = [1 0; 0 1];

% Assuming process noise is white for each component, covariance matrix of
% process noise must be a diagonal matrix, with noise variance for each
% component in each position of the diagonal.
Q1 = beta1*[gyroVar 0; 0 accVar];
Q2 = beta2*[gyroVar 0; 0 accVar];
% Q = [gyroVar 0; 0 accVar];

% NEED TO JUSTIFY THE FORM OF MATRIX A.
A = [1 -dt; 0 1];

% Bias is not needed to compute angle estimate.
C = [1 0];

% Initializing state vector
X = [Ini; gyro(1)];
x1 = zeros(len, 1);
x2 = zeros(len, 1);

for i=1:1:len

        if marker(i)==0
            Rk=Rk1;
            Q=Q1;
        end
        if marker(i)==1
            Rk=Rk2;
            Q=Q2;
        end
        % Predition
        X(1) = X(1) + (gyro(i)-X(2))*dt;
        X(2) = X(2);
        P = A*P*A' + Q;

        % Updating estimation
        v = acc(i) - C*X;
        Sk = C*P*C' + Rk;
        K = (P*C')/Sk;
        X = X + K*v;
        P = P - K*Sk*K';
        x1(i) = X(1);
        x2(i) = X(2);
        P_output(:,:,i)=P;
end
end