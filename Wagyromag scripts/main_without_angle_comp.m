% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% **************************** WAGYROMAG **********************************
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% 
% This file includes a complete toolkit to load, calibrate and process the
% data gathered with the Wagyromag MIMU.
%
% -------------------------------------------------------------------------
% * Author:   Alberto Olivares Vicente.
% * Entity:   Universidad de Granada. 
% * Version:  1.0.
% * Last modification: 20/01/2014.
% -------------------------------------------------------------------------


clear all, close all;

% -------------------------------------------------------------------------
% 0) GENERAL CONFIGURATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -------------------------------------------------------------------------
% Load Wagyromag's library.
wag = wagLibrary;

% Set Wagyromag's ID:
wag_id = 15;

% Set sampling frequency.
f = 40;               

% Set value of the magnitude of the gravity vector in the location in which
% data were gathered. (In our case: Granada, Spain, 37?10'4''N 3?36'3''O, 
% 738 meters over sea level). 
g = 9.797024;

% -------------------------------------------------------------------------
% 1) LOAD AND CALIBRATE RAW DATA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -------------------------------------------------------------------------
% This section reads the following data from  the .CSV file. Raw data are 
% stored in the following vectors of size (1 x N), where N is the total
% number of samples. 
%  - ax: Acceleration measured along X axis. Raw values are within the 
%        [0-1023] range. 
%  - ay: Acceleration measured along Y axis. Raw values are within the 
%        [0-1023] range.
%  - az: Acceleration measured along Z axis. Raw values are within the 
%        [0-1023] range.
%  - gx: Angular rate measured along X axis. Raw values are within the 
%        [0-1023] range.
%  - gy: Angular rate measured along Y axis. Raw values are within the 
%        [0-1023] range.
%  - gz: Angular rate measured along Z axis. Raw values are within the 
%        [0-1023] range.
%  - hx: Magnetic field measured along X axis. Raw values are within the 
%        [0-1023] range.
%  - hy: Magnetic field measured along Y axis. Raw values are within the 
%        [0-1023] range.
%  - hz: Magnetic field measured along Z axis. Raw values are within the 
%        [0-1023] range.
%  - Angle: Actual orientation. This vector is used as a reference to
%           check the accuracy of the different methods.

% -------------------------------------------------------------------------
% 1.1) Browse for data file and load data.
% -------------------------------------------------------------------------
[filename, pathname] = uigetfile('*.CSV', 'Select data file');
completepath = fullfile(pathname, filename);
[ax, ay, az, gx, gy, gz, hx, hy, hz] = wag.load_data(completepath); 

% -------------------------------------------------------------------------
% 1.2) Calibrate raw data. 
% -------------------------------------------------------------------------
% Once data are loaded, we proceed to calibrate them, i.e. we transform the
% raw units [0-1023] into meaningful physical units (g, deg/s and Gauss). 

% Calibrate ACCELERATION.
[axC, ayC, azC] = wag.calibrate(ax, ay, az, wag_id, 'acc'); 

% Truncate noise (hardcore style).
axC(abs(axC(:))<1e-6) = 0;      
ayC(abs(ayC(:))<1e-6) = 0;      
azC(abs(azC(:))<1e-6) = 0;  

% Transform units (from 'm/s^2' to 'g').
axC = axC / g;                    
ayC = ayC / g;                   
azC = azC / g;                   

% Axis adjustment (fitting the MIMU's body frame to the inertial frame). In
% this case axes X and Z are interchanged. 
axCaux = azC;       
azCaux = axC;       
axC = axCaux;         
azC = azCaux;

% Calibrate MAGNETIC FIELD. 
[hxC, hyC ,hzC] = wag.calibrate(hx, hy, hz, wag_id, 'mag'); 

% Axis adjustment (fitting the MIMU's body frame to the inertial frame). In
% this case, the sense of axes Y and Z is inverted. 
hyC = -hyC;                       
hzC = -hzC;                       

% Calibrate ANGULAR RATE. 
[gxC, gyC, gzC] = wag.calibrate(gx, gy, gz, wag_id, 'gyro'); 

% Compensate remaining deviation (hardcore style).
gxC = gxC - gxC(1);          
gyC = gyC - gyC(1);         
gzC = gzC - gzC(1);  

% Transform units from deg/s to rad/s
gxC =  gxC * pi/180;          
gyC =  gyC * pi/180;         
gzC = -gzC * pi/180;         

% Build the acceleration, angular rate and magnetic field matrices
Accelerometer = [axC ayC azC]; 
Gyroscope = [gxC gyC -gzC];    
Magnetometer = [hxC' hyC' hzC']; 

% Build the time instants vector.
time = zeros(1, length(axC))';  
for j = 2:length(axC)
    time(j) = time(j-1) + 1 / f;
end