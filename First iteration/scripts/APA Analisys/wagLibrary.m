function wag = wagLibrary

% FUNCTION WAGLIBRARY contains the references to all necesarry functions to 
% load, extract, calibrate, process and plot data gathered using Wagyromag.
% Additionally, it contains all the necessary functions for the study of
% gating in orientation estimation algorithms.

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

% Data loading and calibration functions.
wag.load_data = @load_data;
wag.calibDigAng = @calibDigAng;
wag.calibrateAcc = @calibrateAcc;
wag.calibrateMag = @calibrateMag;
wag.calibrateGyro = @calibrateGyro;

% Functions to estimate orientation using the acceleration and magnetic
% field decomposition.
wag.pitch_roll_decomp = @pitch_roll_decomp;
wag.yaw_decomp = @yaw_decomp;

% QUEST algorithm functions.
wag.quest = @quest;
wag.optimizeQUEST = @optimizeQUEST;

% Regular Kalman filter functions.
wag.optimizeKF = @optimizeKF;
wag.fusionKF = @fusionKF;

% Extended Kalman filter functions.
wag.optimizeEKF = @optimizeEKF;
wag.fusionEKF = @fusionEKF;

% Extended Kalman filter + QUEST functions.
wag.optimizeEKF_QUEST = @optimizeEKF_QUEST;
wag.fusionEKF_QUEST = @fusionEKF_QUEST;

% Gated Kalman filter functions.
wag.optimizeGKF = @optimizeGKF;
wag.fusionGKF = @fusionGKF;

% Gated Extended Kalman filter functions.
wag.optimizeGEKF = @optimizeGEKF;
wag.fusionGEKF = @fusionGEKF;

% Gated Extended Kalman filter + QUEST functions.
wag.optimizeGEKF_QUEST = @optimizeGEKF_QUEST;
wag.fusionGEKF_QUEST = @fusionGEKF_QUEST;

% Auxiliary functions.
wag.quat_to_euler = @ quat_to_euler;
wag.create_figures = @ create_figures;
wag.compute_rmse = @compute_rmse;
wag.build_int_markers = @build_int_markers;
wag.improvement_perc = @improvement_perc;

% Activity detectors.
wag.fsd = @fsd;
wag.ltsd = @ltsd;
wag.compEstMark = @compEstMark;

end
% END OF WAGLIBRARY FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [ax, ay, az, gx, gy, gz, hx, hy, hz, temp, bat, angle] = ...
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
% - INPUT PARAMETERS:
%   |_'pathname' (string): Path to the folder where the data file is
%     stored.
%   |_'filename' (string): Name of the data file.
% 
% - OUTPUT PARAMETERS:
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

% Read the complete data matrix (all columns)
alldata = dlmread(completepath, ';', 3, 1); 

% Extract data columns.
ax = alldata(:, 1);
ay = alldata(:, 2);
az = alldata(:, 3);
gx = alldata(:, 4);
gy = alldata(:, 5);
gz = alldata(:, 6);
hx = alldata(:, 7);
hy = alldata(:, 8);
hz = alldata(:, 9);
temp = alldata(:, 10);
bat = alldata(:, 11);
angle = alldata(:, 12);

end
% END OF LOAD_DATA FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


function true_ang = calibDigAng(dig_ang)

% FUNCTION CALIBDIGANG transforms raw reference angle data into degrees
% using a look up table that should be already stored in the system.
%
% - INPUT PARAMETERS:
%   |_'dig_ang' (vector): Raw reference angle.
%
% - OUTPUT PARAMETERS:
%   |_'true_ang' (vector): Reference angle in degrees.
% -------------------------------------------------------------------------
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
% Load look up table.
load('data/LookUpTableSplineMapRaw2Deg.mat');

% Map raw data to calibrated data.
true_ang = zeros(1, length(dig_ang));
for i = 1:length(dig_ang)
    if dig_ang(i) == 0
        true_ang(i) = 0;
    else 
        true_ang(i) = y_spline(dig_ang(i));
    end
end

end
% END OF CALIBDIGANG FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [axC, ayC, azC] = calibrateAcc(ax, ay, az, Ka, Ra, ba)

% FUNCTION CALIBRATEACC applies the calibration parameters (Ka, Ra and ba) 
% to the raw accelerometer signals. 
%
% - INPUT PARAMETERS:
%   |_'ax' (vector): Raw acceleration along X-axis.
%   |_'ay' (vector): Raw acceleration along Y-axis.
%   |_'az' (vector): Raw acceleration along Z-axis.
%   |_'Ka' (3x3 matrix): Scale factor matrix.
%   |_'Ra' (3x3 matrix): Misalignment/orientation matrix.
%   |_'ba' (1x3 vector): Accelerometer biases.
%
% - OUTPUT PARAMETERS:
%   |_'axC' (vector): Calibrated acceleration along X-axis.
%   |_'ayC' (vector): Calibrated acceleration along Y-axis.
%   |_'azC' (vector): Calibrated acceleration along Z-axis. 
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
% Check input parameters.
len = length(ax);
if length(ay) ~= len || length(az) ~= len
    error(['Input parameters ''accX'', ''accY'' and ''accZ'' must be ' ...
        'the same length.']);
end

% Apply the calibration formula which contains the calibration parameters.
aC = [ones(1, length(ax))' ones(1, length(ay))' ones(1, length(az))'];

for i = 1:length(ax)
    a = [ax(i) ay(i) az(i)]';
    
    % Calibration formula.
    aC(i, :) = (inv(Ra) * inv(Ka) * (a - ba'))';
end

% Extract vectors from matrix.
axC = aC(:, 1);  
ayC = aC(:, 2);  
azC = aC(:, 3);

end
% END OF CALIBRATEACC FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [hxC, hyC, hzC] = calibrateMag(hx, hy, hz, T, S, b)

% FUNCTION CALIBRATEMAG applies the calibration parameters (T, S and b) to 
% the raw magnetic field signals. 
% 
% - INPUT PARAMETERS:
%   |_'hx' (vector): Raw magnetic field along X-axis.
%   |_'hy' (vector): Raw magnetic field along Y-axis.
%   |_'hz' (vector): Raw magnetic field along Z-axis.
%   |_'T' (3x3 matrix): Misalignment matrix.
%   |_'S' (3x3 matrix): Scale factor matrix.
%   |_'b' (3x1 vector): bias.
%
% - OUTPUT PARAMETERS:
%   |_'hxC': Calibrated magnetic field along X-axis.
%   |_'hyC': Calibrated magnetic field along Y-axis.
%   |_'hzC': Calibrated magnetic field along Z-axis. 
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%   |_ Gonzalo Ruiz-García:
%       * Entity:   Department of Computer Architecture and Computer
%                   Technology, University of Granada, Spain.
%       * Contact:  gruiz@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares & Gonzalo Ruiz-García. 
%  All rights reserved.
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

% Check input parameters.
len = length(hx);
if (length(hy) ~=len || length(hz) ~= len)
    error(['Input parameters ''hx'', ''hy'' and ''hz'' must ' ...
        'be the same length.']);
end
if size(T, 1) ~= 3 || size(T, 2) ~= 3
    error('Input parameter ''T'' must be a 3x3 matrix.');
end
if size(S, 1) ~= 3 || size(S, 2) ~= 3
    error('Input parameter ''S'' must be a 3x3 matrix.');
end
if size(b, 1) ~= 3 || size(b, 2) ~= 1
    error('Input parameter ''b'' must be a 3x1 matrix.');
end

% Apply the calibration formula which contains the calibration parameters.
hxC = zeros(1, len);
hyC = zeros(1, len);
hzC = zeros(1, len);

for i = 1:1:len   
    % Calibration formula.
    U = inv(T) * inv(S) * ([hx(i); hy(i); hz(i)] - b);
    hxC(i) = U(1);
    hyC(i) = U(2);
    hzC(i) = U(3);
end
end
% END OF CALIBRATEMAG FILE

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [gxC, gyC, gzC] = calibrateGyro(gx, gy, gz, Kg, Rg, bg)

% FUNCTION CALIBRATEGYRO applies the calibration parameters (Kg, Rg and bg) 
% to the raw angular rate signals. 
%
% - INPUT PARAMETERS:
%   |_'gx' (vector): Raw angular rate along X-axis.
%   |_'gy' (vector): Raw angular rate along Y-axis.
%   |_'gz' (vector): Raw angular rate along Z-axis.
%   |_'Kg' (3x3 matrix): Scale factor matrix.
%   |_'Rg' (3x3 matrix): Misalignment/orientation matrix.
%   |_'bg' (3x1 vector): bias.
%
% - OUTPUT PARAMETERS:
%   |_'gxC': Calibrated angular rate along X-axis.
%   |_'gyC': Calibrated angular rate along Y-axis.
%   |_'gzC': Calibrated angular rate along Z-axis. 
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%   |_ Gonzalo Ruiz-García:
%       * Entity:   Department of Computer Architecture and Computer
%                   Technology, University of Granada, Spain.
%       * Contact:  gruiz@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares & Gonzalo Ruiz-García. 
%  All rights reserved.
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

% Check input parameters.
len = length(gx);
if (length(gy) ~= len || length(gz) ~= len)
    error(['Input parameters ''gx'', ''gy'' and ''gz'' must be the ' ...
        'same length.']);
end
if size(Kg, 1) ~= 3 || size(Kg, 2) ~= 3
    error('Input parameter ''Kg'' must be a 3x3 matrix.');
end
if size(Rg, 1) ~= 3 || size(Rg, 2) ~= 3
    error('Input parameter ''Rg'' must be a 3x3 matrix.');
end
if size(bg, 1) ~= 1 || size(bg, 2) ~= 3 
    error('Input parameter ''bg'' must be a 1x3 vector.');
end

% Apply the calibration formula which contains the calibration parameters.
gC = [ones(1, len)' ones(1, len)' ones(1, len)'];      
for i = 1:len
    g = [gx(i) gy(i) gz(i)]' - bg';
    gC(i, :) = (inv(Rg) * inv(Kg) * (g))';
end

% Extract vectors from matrix.
gxC = gC(:, 1);
gyC = gC(:, 2);
gzC = gC(:, 3);

end
% END OF CALIBRATEGYRO FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [pitch, roll] = pitch_roll_decomp(ax, ay, az, units)

% FUNCTION PITCH_ROLL_ACC Computes pitch and roll angles using the
% decomposition of the acceleration measured with a triaxial accelerometer. 
%
% - INPUT PARAMETERS:
%    |_'ax': vector containing the linear accelerations along the cartesian
%          X-axis.
%    |_'ay': vector containing the linear accelerations along the cartesian
%          Y-axis.
%    |_'az': vector containing the linear accelerations along the cartesian
%          Z-axis.
%    |_'units': string containing the units of angles to return. Set 'deg' 
%               or 'rad'.
%
% - OUTPUT PARAMETERS:
%    |_'pitch': vector containing the pitch calculated angle values given 
%               the 'ay' and 'az' components of linear accelerations. Pitch
%               is assumed to be the rotation angle about the Y-axis.
%    |_'roll':  vector containing the roll calculated angle values given 
%               the 'ax' and 'az' components of linear accelerations. Roll
%               is asumed to be the rotation angle about the X-axis.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%   |_ Gonzalo Ruiz-García:
%       * Entity:   Department of Computer Architecture and Computer
%                   Technology, University of Granada, Spain.
%       * Contact:  gruiz@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares & Gonzalo Ruiz-García. 
%  All rights reserved.
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

% Check input parameters.
if (length(ax) ~= length(ay) || length(ay) ~= length(az) || ...
    length(az) ~= length(ax))
    error('Input vectors ax, ay and az must be the same length.');
end

if ~(strcmp(units, 'deg') || strcmp(units, 'rad'))
    error('Specified units are not correct. Set ''deg'' or ''rad''.');
end

% Compute pitch and roll. Pitch and roll values are obtained (in rad). Both
% pitch and roll values have different expressions depending on the 
% quadrant they are placed.
pitch = atan2(sqrt(ay .^ 2 + az .^ 2), -ax);
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
% - INPUT PARAMETERS: 
%    |_ 'pitch': vector containing the a priori known pitch angle values.
%    |_ 'roll': vector containing the a priori known roll angle values.
%    |_ 'magX': vector containing the a priori known magnetic field 
%               component along X-axis.
%    |_ 'magY': vector containing the a priori known magnetic field 
%               component along Y-axis.
%    |_ 'magZ': vector containing the a priori known magnetic field 
%               component along Z-axis.
%    |_ 'plotGraphics': string containing 'yes' or 'no' to specify if 
%                       results are to be plotted.
%
% - OUTPUT PARAMETERS:
%    |_ 'yaw': vector containing the calculated yaw angles, considering 
%              both local magnetic field and pitch and roll values.
%    |_ 'magXh': vector containing the projections of magX over the XY 
%                cartesian plane. These values are used to compute YAW 
%                output values.
%    |_ 'magYh': vector containing the projections of magY over the XY 
%                cartesian plane. These values are used to compute YAW 
%                output values.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%   |_ Gonzalo Ruiz-García:
%       * Entity:   Department of Computer Architecture and Computer
%                   Technology, University of Granada, Spain.
%       * Contact:  gruiz@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares & Gonzalo Ruiz-García. 
%  All rights reserved.
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

% Check input parameters.
len = length(pitch);
if (length(roll) ~= len || length(magX) ~= len || length(magY) ~= ...
        len || length(magZ) ~= len)
    error('All input arguments must be the same length.');
end

% Compute Yaw Angle. Calculation of horizontal magnetic X coordinate 
% component (magXh) and horizontal magnetic Y coordinate component (magYh).
% This is done by derotating the measured magnetic field in axes X and Y by
% the pitch and roll angles of the body. That is, we find the projection of
% the magnetic field in axes X and Y in the XY plane. 
magXh = magX .* cos(pitch) + magY .* sin(roll) .* sin(pitch) - ...
        magZ .* cos(roll) .* sin(pitch);
magYh = magY .* cos(roll) + magZ .* sin(roll);
yaw = atan((-magYh) ./ magXh);

% Quadrant compensations. This part may need to be changed depending on the
% body reference frame of the MIMU. 
for i = 1 : length(yaw)
    if magXh(i) > 0 && magYh(i) == 0
        yaw(i) = pi / 2;
    end
    if magXh(i) < 0 && magYh(i) == 0
        yaw(i) = -pi / 2;
    end
    if magXh(i) < 0 && magYh(i) > 0
        yaw(i) = yaw(i) - pi;
    end
    if magXh(i) < 0 && magYh(i) < 0
        yaw(i) = yaw(i) + pi;
    end
    
    while yaw(i) > pi 
        yaw(i) = yaw(i) - 2 * pi; 
    end
    while yaw(i) < -pi 
        yaw(i) = yaw(i) + 2 * pi; 
    end
end

end
% END OF YAW_DECOMP FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [roll_QUEST, pitch_QUEST, yaw_QUEST, q_quest] = quest(axC, ayC,...
    azC, hxC, hyC, hzC, a1, a2)

% FUNCTION QUEST carries out the QUEST method proposed by by M. D. Shuster 
% and S. D. Oh in the work entitled: "Three-axis attitude determination 
% from vector observations". Journal of Guidance Control and Dynamics, 
% 4(1):70–77, 1981.
%
% - INPUT PARAMETERS:
%   |_ 'axC': vector containing acceleration gathered along X-axis.
%   |_ 'ayC': vector containing acceleration gathered along Y-axis.
%   |_ 'azC': vector containing acceleration gathered along Z-axis.
%   |_ 'hxC': vector containing magnetic component along X-axis.
%   |_ 'hyC': vector containing magnetic component along Y-axis.
%   |_ 'hzC': vector containing magnetic component along Z-axis.
%   |_ 'a1': Weight of the measured acceleration in the computation.
%   |_ 'a2': Weight of the measured magnetic field in the computation.    
%
% - OUTPUT PARAMETERS:
%   |_ 'roll_QUEST': Computed roll.
%   |_ 'pitch_QUEST': Computed pitch.
%   |_ 'yaw_QUEST': Computed yaw.
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

% Set the gravity vector in earth frame according to the predefined
% coordinate system.
v1 = [-1 0 0]; 

% Set magnetic field vector in Earth Frame at the location of the Lab 
% (Granada, Spain).
v2 = [0.27195  -0.00634  -0.33741]; 
a = 1 / (norm(cross(v1, v2), 2) ^ 2);
b = (v1 .* v2) / (norm(cross(v1, v2), 2) ^ 2);

for i = 1:length(axC)
    w1 = [axC(i) ayC(i) azC(i)];
    w2 = [hxC(i) hyC(i) hzC(i)];
    Y = 1 / 2 * (a1 * cross(w1, v1) + a2 * cross(w2, v2) + a1 * ((b .* ...
        w1 - a .* w2) .* cross(v1, v2)) .* v1 + a2 * ((a .* w1 - ...
        b .* w2) .* cross(v1, v2)) .* v2);
    q_quest(i, :) = (1 / sqrt(1 + norm(Y, 2) ^ 2) * [Y'; 1])';
    clear eulerAngles;
    eulerAngles = quat_to_euler(q_quest(i, :));
    roll_QUEST(i) = eulerAngles(1);
    pitch_QUEST(i) = -eulerAngles(2);
    yaw_QUEST(i) = eulerAngles(3);
end

end
% END OF QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeQUEST(ax, ay, az, hx, hy, hz, ...
    ref_angle, rotation_angle, p0_QUEST, rmse_off)

% FUNCTION OPTIMIZEQUEST uses the ANMS Algorithm to optimize the parameters
% of the QUEST ALGORITHM by minimizing the error between the actual 
% orientation angle and the one estimated with the QUEST algorithm. 
%
% - INPUT VARIABLES:
%   |_ 'ax': Acceleration measured along X axis.
%   |_ 'ay': Acceleration measured along Y axis.
%   |_ 'az': Acceleration measured along Z axis.
%   |_ 'hx': Magnetic field measured along X axis.
%   |_ 'hy': Magnetic field measured along Y axis.
%   |_ 'hz': Magnetic field measured along Z axis.
%   |_ 'ref_angle': Reference orientation angle (actual value).
%   |_ 'rotation_angle': Rotation angle around which data were gathered
%                        (possible values are 'pitch, 'roll' or 'yaw');
%   |_ 'p0_QUEST': Vector containing initial guess of the parameters.
%
% - OUTPUT VARIABLES:
%   |_ 'xmin': Value of the parameters minimizing the error function (found
%              by the ANMS algorithm).
%   |_ 'fmin': Minimum value of the error function.
%   |_ 'ct':   Number of ANMS iterations until the minimum was found.
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
% 1) Set global variables.
% -------------------------------------------------------------------------
% The variables must be declared as global in order to be used inside the
% eofQUEST function.
global acc_x;           acc_x = ax;
global acc_y;           acc_y = ay;
global acc_z;           acc_z = az;
global mag_x;           mag_x = hx;
global mag_y;           mag_y = hy;
global mag_z;           mag_z = hz;
global true_angle;      true_angle = ref_angle;
global chosen_angle;    chosen_angle = rotation_angle;
global rmse_offset;     rmse_offset = rmse_off;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp('Optimizing parameters of QUEST (this may take a while) ...');

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-4;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin, fmin, ct] = ANMS(@eofQUEST, p0_QUEST, tol, max_feval);

end

% END OF OPTIMIZEQUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofQUEST(p)

% FUNCTION EOFQUEST Computes the error function to be minimized: RMSE
% between the actual angle and the estimated orientation angle computed
% with the QUEST Algorithm. 
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
%
% 1) Set variables.
% -------------------------------------------------------------------------
% Define global variables.
global acc_x;
global acc_y;
global acc_z;
global mag_x;
global mag_y;
global mag_z;
global true_angle;
global chosen_angle;
global rmse_offset;

% Get variables from input vector.
a1 = p(1);
a2 = p(2);

% Call the QUEST algorithm to estimate the orientation angle.
[roll_QUEST, pitch_QUEST, yaw_QUEST] = quest(acc_x, acc_y, acc_z, ...
mag_x, mag_y, mag_z, a1, a2);

% Select the chosen angle.
if strcmp(chosen_angle,'pitch')
    angle_QUEST = pitch_QUEST;
elseif strcmp(chosen_angle,'roll')
    angle_QUEST = roll_QUEST;
elseif strcmp(chosen_angle,'yaw')
    angle_QUEST = yaw_QUEST;
end

% Define the error function to be minimized. 
F = compute_rmse(true_angle, 180 / pi* angle_QUEST', rmse_offset);

end
% END OF EOFQUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function euler = quat_to_euler(quatInput)

% FUNCTION QUAT_TO_EULER Extracts the euler angles from a given orientation
% quaternion. 
%
% - INPUT PARAMETERS:
%    |_ 'quatInput': Nx4 matrix containing a quaternion per row. Each row 
%                    will be converted to a 1x3 matrix containing Euler's 
%                    corresponding angles.
%
% - OUTPUT PARAMETERS:
%    |_ 'euler': Nx3 matrix containing corresponding Euler's angles for 
%                each input quaternion. Angles are as follows:
%                   - euler(:,1): rotation about the X-axis.
%                   - euler(:,2): rotation about the Y-axis.
%                   - euler(:,3): rotation about the Z-axis.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%   |_ Gonzalo Ruiz-García:
%       * Entity:   Department of Computer Architecture and Computer
%                   Technology, University of Granada, Spain.
%       * Contact:  gruiz@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/21/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Alberto Olivares & Gonzalo Ruiz-García. 
%  All rights reserved.
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

% 0) Check input parameters.
% -------------------------------------------------------------------------
if size(quatInput, 2) ~= 4
    error('Input parameter ''quatInput'' must be a Nx4 matrix.');
end

% 1) Set variables.
% -------------------------------------------------------------------------
len = size(quatInput, 1);
euler = zeros(len, 3);

% 2) Apply transformations.
% -------------------------------------------------------------------------
for i = 1 : len
    q0 = quatInput(i, 1);
    q1 = quatInput(i, 2);
    q2 = quatInput(i, 3);
    q3 = quatInput(i, 4);
    euler(i, 1) = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + ...
        q2 * q2));
    euler(i, 2) = asin(2 * (q0 * q2 - q3 * q1));
    euler(i, 3) = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + ...
        q3 * q3));
end
end
% END OF QUAT_TO_EULER FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function create_figures(order, method_name, x_axis, y_axis, format)

% FUNCTION CREATE_FIGURES Plots the orientation estimates. The function
% updates the legend automatically when a new signal is added to the plot.
% There is only a distinction between the first time the function is called
% and the subsequent executions. 
% 
% - INPUT PARAMETERS:
%   |_ 'order': Order in which the function is called. 
%               Possible values:
%                   - 'first': To indicate that the function is called for
%                   the first time.
%                   - 'other': Any other time but the first.
%   |_ 'method_name': Name of the orientation estimation algorithm which
%                     will be displayed in the legend.
%   |_ 'x_axis': Signal to be displayed in the X axis. Usually the time
%                signal.
%   |_ 'y_axis': Signal to be displayed in the Y axis. Usually the angle
%                signal.
%
% - OUTPUT PARAMETERS: This function has no output parameters.
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

% If the function is called for the first time, plot the signals normally.
if strcmp(order, 'first')
    plot(x_axis, y_axis, format, 'LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    legend(method_name);
end

% If the function has already been called at least once, plot the signals,
% get the names in the legend and update them. 
if strcmp(order, 'other')
    p2 = plot(x_axis, y_axis, format, 'LineWidth', 2);
    hold on
    ch = get(gca, 'children');
    ls = get(legend, 'String');
    if isempty(ls)
        legend(p2, method_name);
        xlabel('Time(s)')
        ylabel('Angle (deg)')
    else
        legend([p2 ch(1)], {char(ls), method_name});
        legend('-DynamicLegend');
        xlabel('Time(s)')
        ylabel('Angle (deg)')
    end 
end
end
% END OF CREATE_FIGURES FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function rmse = compute_rmse(reference, estimation, rmse_offset)

% FUNCTION COMPUTE_RMSE Computes the RMSE (Root Mean Squared Error)of the 
% estimated orientation angle with respect to the ground truth (reference).
% The rmse is computed from the Xth sample to the Nth sample, where X is an 
% initial offset and N is the total length of the signal. This is done to 
% allow the adaptive algorithms to converge. 
% 
% - INPUT PARAMETERS:
%   |_ 'reference': Reference signal acting as the ground truth.
%   |_ 'estimation': Signal containing the orientation angle estimate.
%   |_ 'rmse_offset': Initial offset of the RMSE computation.
%
% - OUTPUT PARAMETERS:
%   |_ 'rmse': Computed RMSE.
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

rmse = sqrt(mean((reference(rmse_offset : end) -  ...
    estimation(rmse_offset : end)') .^ 2));
        
end
% END OF CREATE_FIGURES FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeKF(observation, ang_vel, frequency, ...
    ref_angle, p0, rmse_off)

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
global gyro;        gyro = ang_vel;
global obs;         obs = observation;
global true_angle;  true_angle = ref_angle;
global freq;        freq = frequency;
global rmse_offset; rmse_offset = rmse_off;

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
global gyro;
global obs;
global freq;
global true_angle;
global rmse_offset;

% 2) Extract the first parameter to be minimized (measurement noise 
% variance gain). 
alpha = p(1); 
beta = p(2);

% 3) Estimate the orientation angle using the Kalman Filter.
angle_KF = fusionKF(gyro, obs, freq, alpha, beta);

% 4) Compute the error function.
F = compute_rmse(true_angle, 180 / pi * angle_KF, rmse_offset);
end
% END OF EOFKALMAN FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [x1, x2, v_output, K_output] = fusionKF(gyro, obs, freq, ...
    alpha, beta)

% FUNCTION FUSIONKF Applies a Kalman Filter sensor fusion approach to
% estimate the orientation of a body using the acceleration, magnetic field
% and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'gyro': vector containing the gyroscope signal for a determined 
%               axis. 
%   |_ 'obs': vector containing the angle values obtained via the linear
%             acceleration relations.
%   |_ 'freq': sampling frecuency. Must be real positive.
%   |_ 'alpha': Weighing factor which multiplies the variance of the
%               measurement noise. It is a parameter which tunes the
%               filter.
%   |_ 'beta': Weighing factor which multiplies the variance of the process
%              noise. It is a parameter which tunes the filter.
%   
% - OUTPUT PARAMETERS:
%   |_ 'x1': estimated orientation angle. (First element of the state
%            vector).
%   |_ 'x2': estimated gyroscope bias. (Second element of the state
%            vector).
%   |_ 'v_output': Vector containing the state vector corrected 'a
%                  posteriori' by the observation.
%   |_ 'K_output': Vector containing the value of the Kalman gain for each
%                  time instant.
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

% 1) Definition of variables.
% -------------------------------------------------------------------------
obsVar = var(obs);
gyroVar = var(gyro);
len = length(gyro);

%Sampling period.
dt = 1/freq;

% Initialization of covariance matrix.
P = [1 0; 0 1];

% Set measurement variance.
Rk = alpha * obsVar;

% Assuming process noise is white for each component, covariance matrix of
% process noise must be a diagonal matrix, with noise variance for each
% component in each position of the diagonal.
Q = [beta * gyroVar 0; 0 beta * obsVar];

% Set state transition matrix.
A = [1 -dt; 0 1];

% Set measurement matrix. 
C = [1 0];

% Initialize state vector.
X = [0; 0];
x1 = zeros(len, 1);
x2 = zeros(len, 1);

% 2) Body of the Kalman Filter.
% -------------------------------------------------------------------------
for i = 1 : len
    
    % Prediction phase.
    X(1) = X(1) + (gyro(i) - X(2)) * dt;
    X(2) = X(2);
    P = A * P * A' + Q;

    % Update phase.
    v = obs(i) - C * X;
    v_output(i) = v;
    Sk = C * P * C' + Rk;
    K = (P * C') / Sk;
    K_output(i, :) = K';
    X = X + K * v;
    P = P - K * Sk * K';
    x1(i) = X(1);
    x2(i) = X(2);   
end

end
% END OF FUSIONKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeEKF(axC, ayC, azC, gxC, gyC, gzC, ...
    frequency, ref_angle, p0, rmse_off, selected_angle)

% FUNCTION OPTIMIZEEKF uses the ANMS algorithm to find the optimal 
% parameters of the Extended Kalman Filter which minimize the error between
% the actual orientation angle and the one estimated with EKF. 
% 
% - INPUT PARAMETERS:
%   |_ 'axC': Calibrated acceleration measured along the X-axis.
%   |_ 'ayC': Calibrated acceleration measured along the Y-axis.
%   |_ 'azC': Calibrated acceleration measured along the Z-axis.
%   |_ 'gxC': Calibrated angular rate measured along the X-axis.
%   |_ 'gyC': Calibrated angular rate measured along the Y-axis.
%   |_ 'gzC': Calibrated angular rate measured along the Z-axis.
%   |_ 'frequency': Sampling frequency.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'rmse_off': Offset in the RMSE computation.
%   |_ 'selected_angle': Orientation angle being estimated. Possible values
%                        are 'pitch', 'roll' and 'yaw'.
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

% 1) Set variables.
% -------------------------------------------------------------------------
global gyro;        gyro = [gxC gyC gzC];
global acc;         acc = [azC ayC axC];
global freq;        freq = frequency;
global true_angle;  true_angle = ref_angle;
global rmse_offset; rmse_offset = rmse_off;
global chosen_angle;chosen_angle = selected_angle;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp(['Optimizing parameters of Extended Kalman Filter (it may take a '...
     'while)...']);

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin,fmin,ct] = ANMS(@eofEKF, p0, tol, max_feval);
 
end
% END OF OPTIMIZEEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofEKF(p)

% FUNCTION EOFEKF Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the Extended Kalman Filter. 
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
global gyro;
global acc;
global true_angle;
global freq;
global chosen_angle;
global rmse_offset;

% Extract the value of the first parameter (variance of gyroscope bias).
sigma_w = p(1);

% Extract the value of the second parameter (variance of the observation).
sigma_obs = p(2);

% Set sampling period. 
dt = 1 /  freq;

% Initialize orientation quaternion. 
ini = [1 0 0 0];

% Call the EKF routine to compute the orientation quaternion.
q = fusionEKF(-acc, gyro, dt, sigma_w, sigma_obs, ini);

% Transform the orientation quaternion to Euler angles.
roll_EKF = zeros(1, length(q));
pitch_EKF = zeros(1, length(q));
yaw_EKF = zeros(1, length(q));
for i = 1 : length(q)
     eulerAngles = quat_to_euler(q(i, :));
     roll_EKF(i) = eulerAngles(1);
     pitch_EKF(i) = eulerAngles(2);
     yaw_EKF(i) = eulerAngles(3);
end

% Select the chosen angle. 
if strcmp(chosen_angle, 'pitch')
    angle_EKF = pitch_EKF;
elseif strcmp(chosen_angle, 'roll')
    angle_EKF = roll_EKF;
elseif strcmp(chosen_angle, 'yaw')
    angle_EKF = yaw_EKF;
end

% Define the error function to be minimized.
F = compute_rmse(true_angle, 180 / pi * angle_EKF', rmse_offset);

end
% END OF EOFEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [q, p1, p2, p3, p4, p5, p6, p7, wb] = fusionEKF(a, w, dt, ...
    sigma_w, sigma_obs, quat_ini)

% FUNCTION FUSIONEKF Applies a Kalman Filter sensor fusion approach to
% estimate the orientation of a body using the acceleration, magnetic field
% and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'a': Acceleration matrix.
%   |_ 'w': Angular velocity matrix.
%   |_ 'dt': Sampling period.
%   |_ 'sigma_w': Variance of the process noise. This is a tuning parameter
%                 of the filter.
%   |_ 'sigma_obs': Variance of the measurement noise. This is a tuning
%                   parameter of the filter.
%   |_ 'quat_ini': Quaternion determining the initial orientation.
%
% - OUTPUT PARAMETERS:
%   |_ 'q': Estimated orientation quaternion.
%   |_ 'p1': Element (1,1) of the covariance matrix.
%   |_ 'p2': Element (2,2) of the covariance matrix.
%   |_ 'p3': Element (3,3) of the covariance matrix.
%   |_ 'p4': Element (4,4) of the covariance matrix.
%   |_ 'p5': Element (5,5) of the covariance matrix.
%   |_ 'p6': Element (6,6) of the covariance matrix.
%   |_ 'p7': Element (7,7) of the covariance matrix.
%   |_ 'wb': Vector containing the bias of the angular velocities.
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
%  Copyright (c) 2015, Alberto Olivares. All rights reserved. Based on the 
%  code in: http://goo.gl/Yww9Bb
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
% Initialize orientation quaternion.
q = zeros(length(a(:, 1)), 4);

% Initialize gyroscope bias vector. 
wb = zeros(length(a(:, 1)), 3);

% Initialize the process noise covariance matrix.
Q = [0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, sigma_w,   0,          0;
     0, 0, 0, 0, 0,         sigma_w,    0;
     0, 0, 0, 0, 0,         0,          sigma_w];

% Initialize the measurement noise covariance matrix.
R = [sigma_obs  0           0  ;
     0          sigma_obs   0  ;
     0          0           sigma_obs;];

% Initialize the system covariance matrix.
P = eye(length(Q)) * 10000; 

% Initialize the state vector
x = [quat_ini 0 0 0]';

% 2) Body of the EKF.
% -------------------------------------------------------------------------
for i = 1 : length(a(:, 1))
    
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    wxb = x(5);
    wyb = x(6);
    wzb = x(7);
    wx = w(i, 1);
    wy = w(i, 2);
    wz = w(i, 3);
    clear z
    z(1) = a(i, 1);
    z(2) = a(i, 2);
    z(3) = a(i, 3);
    z = z';

    % Populate F jacobian
    F = [              1,  (dt/2)*(wx-wxb),  (dt/2)*(wy-wyb), ...
         (dt/2)*(wz-wzb),       -(dt/2)*q1,       -(dt/2)*q2, -(dt/2)*q3;
        -(dt/2)*(wx-wxb),                1,  (dt/2)*(wz-wzb), ...
        -(dt/2)*(wy-wyb),        (dt/2)*q0,        (dt/2)*q3, -(dt/2)*q2;
        -(dt/2)*(wy-wyb), -(dt/2)*(wz-wzb),                1, ...
         (dt/2)*(wx-wxb),       -(dt/2)*q3,        (dt/2)*q0, (dt/2)*q1;
        -(dt/2)*(wz-wzb),  (dt/2)*(wy-wyb), -(dt/2)*(wx-wxb), ...  
                       1,        (dt/2)*q2,       -(dt/2)*q1, (dt/2)*q0;
                       0,                0,                0, ...
                       0,                1,                0,          0;
                       0,                0,                0, ...
                       0,                0,                1,          0;
                       0,                0,                0, ...
                       0,                0,                0,          1;];
                   
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Predicted state estimate, x = f(x,u).
    x = [q0 - (dt/2) * (-q1*(wx-wxb) - q2*(wy-wyb) - q3*(wz-wzb));
         q1 - (dt/2) * ( q0*(wx-wxb) + q3*(wy-wyb) - q2*(wz-wzb));
         q2 - (dt/2) * (-q3*(wx-wxb) + q0*(wy-wyb) + q1*(wz-wzb));
         q3 - (dt/2) * ( q2*(wx-wxb) - q1*(wy-wyb) + q0*(wz-wzb));
         wxb;
         wyb;
         wzb;];

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);

    % Predicted covariance estimate.
    P = F * P * F' + Q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Normalize accelerometer and magnetometer measurements.
    acc_norm = sqrt(z(1) ^ 2 + z(2) ^ 2 + z(3) ^ 2);
    z(1) = z(1) / acc_norm;
    z(2) = z(2) / acc_norm;
    z(3) = z(3) / acc_norm;

    h = [-2 * (q1 * q3 - q0 * q2);
         -2 * (q2 * q3 + q0 * q1);
         -2 * q0 ^ 2 + 1 - 2 * q3 ^ 2];

    % Measurement residual: y = z - h(x), where h(x) is the matrix that 
    % maps the state onto the measurement.
    y = z - h;

    % The H matrix maps the measurement to the states.
    H = [ 2 * q2, -2 * q3,  2 * q0, -2 * q1, 0, 0, 0;
         -2 * q1, -2 * q0, -2 * q3, -2 * q2, 0, 0, 0;
         -4 * q0,       0,       0, -4 * q3, 0, 0, 0;]; 
     
    % Measurement covariance update.
    S = H * P * H' + R;

    % Calculate Kalman gain.
    K = P * H'/ S;

    % Corrected model prediction.
    x = x + K * y;      

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    % Update state covariance with new knowledge.
    I = eye(length(P));
    P = (I - K * H) * P;  

    % Populate the orientation quaternion matrix.
    q(i,:) = [x(1), x(2), x(3), x(4)];
    
    % Populate the gyroscope bias matrix.
    wb(i,:) = [x(5), x(6), x(7)];
    
    % Get the diagonal elements of the covariance matrix.
    p1(i) = P(1, 1);
    p2(i) = P(2, 2);
    p3(i) = P(3, 3);
    p4(i) = P(4, 4);
    p5(i) = P(5, 5);
    p6(i) = P(6, 6);
    p7(i) = P(7, 7);
end

end
% END OF FUSIONEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeEKF_QUEST(angular_vel, obs, p0, ...
    ref_angle, selected_angle, frequency, rmse_off)

% FUNCTION OPTIMIZEEKF_QUEST uses the ANMS algorithm to find the optimal 
% parameters of the Extended Kalman Filter + QUEST which minimize the error
% between the actual orientation angle and the one estimated with EKF. 
% 
% - INPUT PARAMETERS:
%   |_ 'angular_vel': Angular velocity (measured with the gyroscope).
%   |_ 'obs': Observation of the Kalman Filter.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'selected_angle': Orientation angle being estimated. Possible values
%                        are 'pitch', 'roll' and 'yaw'.
%   |_ 'frequency': Sampling frequency.
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

% 1) Set variables.
% -------------------------------------------------------------------------
global gyro;            gyro = angular_vel;
global quat_obs;        quat_obs = obs;
global true_angle;      true_angle = ref_angle;
global chosen_angle;    chosen_angle = selected_angle;
global freq;            freq = frequency;
global rmse_offset;     rmse_offset = rmse_off;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp(['Optimizing parameters of Extended Kalman Filter (QUEST)',...
    ' (it may take a while) ...']);

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin, fmin, ct] = ANMS(@eofEKF_QUEST, p0, tol, max_feval);

end
% END OF OPTIMIZEEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofEKF_QUEST(p)

% FUNCTION EOFEKF_QUEST Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the EKF-QUEST Algorithm. 

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
global gyro;
global freq;
global quat_obs;
global true_angle;
global chosen_angle;
global rmse_offset;

dt = 1 / freq;
sigma_w = p(1);
sigma_obs = p(2);
quat_ini = quat_obs(1, :);

% 2) Compute the orientation estimation using the EKF-QUEST algorithm.
quat_est = fusionEKF_QUEST(quat_obs, gyro, dt, sigma_w, sigma_obs, ...
    quat_ini);

roll_EKF_QUEST = zeros(1, length(quat_est));
pitch_EKF_QUEST = zeros(1, length(quat_est));
yaw_EKF_QUEST = zeros(1, length(quat_est));
for i = 1 : length(quat_est)
     eulerAngles = quat_to_euler(quat_est(i, :));
     roll_EKF_QUEST(i) = eulerAngles(1);
     pitch_EKF_QUEST(i) = -eulerAngles(2);
     yaw_EKF_QUEST(i) = eulerAngles(3);
end

% 3) Select the chosen angle. 
if strcmp(chosen_angle, 'pitch')
    angle_EKF_QUEST = pitch_EKF_QUEST;
elseif strcmp(chosen_angle, 'roll')
    angle_EKF_QUEST = roll_EKF_QUEST;
elseif strcmp(chosen_angle, 'yaw')
    angle_EKF_QUEST = yaw_EKF_QUEST;
end

% 4) Define the error function to be minimized. 
F = compute_rmse(true_angle, 180 / pi * angle_EKF_QUEST', rmse_offset);

end
% END OF EOFEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [q, wb, y_output] = fusionEKF_QUEST(quat, w, dt, sigma_w, ...
    sigma_obs, quat_ini)

% FUNCTION FUSIONEKF_QUEST Applies a Kalman Filter sensor fusion approach 
% to estimate the orientation of a body using the acceleration, magnetic 
% field and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'quat': Quaternion containing the orientation estimated by the QUEST
%              algorithm. These estimates are used as the observation.
%   |_ 'w': Angular velocity matrix.
%   |_ 'dt': Sampling period.
%   |_ 'sigma_w': Variance of the process noise. This is a tuning parameter
%                 of the filter.
%   |_ 'sigma_obs': Variance of the measurement noise. This is a tuning
%                   parameter of the filter.
%   |_ 'quat_ini': Quaternion determining the initial orientation.
%
% - OUTPUT PARAMETERS:
%   |_ 'q': Estimated orientation quaternion.
%   |_ 'wb': Vector containing the bias of the angular velocities.
%   |_ 'y_output': Measurement residual of every time instant.
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
%  Copyright (c) 2015, Alberto Olivares. All rights reserved. Based on the 
%  code in: http://goo.gl/Yww9Bb
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
% Initialize orientation quaternion.
q = zeros(length(quat(:, 1)), 4);

% Initialize gyroscope bias vector. 
wb = zeros(length(w(:, 1)), 3);

% Initialize the process noise covariance matrix.
Q = [0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, 0,         0,          0;
     0, 0, 0, 0, sigma_w,   0,          0;
     0, 0, 0, 0, 0,         sigma_w,    0;
     0, 0, 0, 0, 0,         0,          sigma_w];

% Initialize the measurement noise covariance matrix.
R = [sigma_obs  0           0          0          0        0        0;
     0          sigma_obs   0          0          0        0        0;
     0          0           sigma_obs  0          0        0        0;
     0          0           0          sigma_obs  0        0        0;
     0          0           0          0          sigma_w  0        0;
     0          0           0          0          0        sigma_w  0
     0          0           0          0          0        0        ...
                                                                  sigma_w];
% Initialize the system covariance matrix.
P = eye(length(Q)) * 10000; 

% Initialize the state vector
x = [quat_ini 0 0 0]';

% 2) Body of the EKF.
% -------------------------------------------------------------------------
for i = 1 : length(quat(:, 1))
    
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    wxb = x(5);
    wyb = x(6);
    wzb = x(7);
    wx = w(i, 1);
    wy = w(i, 2);
    wz = w(i, 3);
    
    % Populate F jacobian
    F = [              1,  (dt/2)*(wx-wxb),  (dt/2)*(wy-wyb), ...
         (dt/2)*(wz-wzb),       -(dt/2)*q1,       -(dt/2)*q2, -(dt/2)*q3;
        -(dt/2)*(wx-wxb),                1,  (dt/2)*(wz-wzb), ...
        -(dt/2)*(wy-wyb),        (dt/2)*q0,        (dt/2)*q3, -(dt/2)*q2;
        -(dt/2)*(wy-wyb), -(dt/2)*(wz-wzb),                1, ...
         (dt/2)*(wx-wxb),       -(dt/2)*q3,        (dt/2)*q0, (dt/2)*q1;
        -(dt/2)*(wz-wzb),  (dt/2)*(wy-wyb), -(dt/2)*(wx-wxb), ...  
                       1,        (dt/2)*q2,       -(dt/2)*q1, (dt/2)*q0;
                       0,                0,                0, ...
                       0,                1,                0,          0;
                       0,                0,                0, ...
                       0,                0,                1,          0;
                       0,                0,                0, ...
                       0,                0,                0,          1;];  
                   
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Predicted state estimate, x = f(x,u).
    x = [q0 - (dt / 2) * (-q1 * (wx - wxb) - q2 * (wy - wyb) - q3 * (wz ...
         - wzb));
         q1 - (dt / 2) * (q0 * (wx - wxb) + q3 * (wy - wyb) - q2 * (wz ...
         - wzb));
         q2 - (dt / 2) * (-q3 * (wx - wxb) + q0 * (wy - wyb) + q1 * (wz ...
         - wzb));
         q3 - (dt / 2) * ( q2 * (wx - wxb) - q1 * (wy - wyb) + q0 * (wz ...
         - wzb));
         wxb;
         wyb;
         wzb;];

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);

    % Predicted covariance estimate.
    P = F * P * F' + Q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h = [q0; q1; q2; q3; 0; 0; 0];

    % Measurement residual: y = z - h(x), where h(x) is the matrix that 
    % maps the state onto the measurement.
    y = [quat(i, :) 0 0 0]' - h;
    y_output(i, :) = y';
    
    % The H matrix maps the measurement to the states.
    H=[1 0 0 0 0 0 0;
       0 1 0 0 0 0 0;
       0 0 1 0 0 0 0;
       0 0 0 1 0 0 0;
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0];
     
    % Measurement covariance update.
    S = H * P * H' + R;

    % Calculate Kalman gain.
    K = P * H' / S;

    % Corrected model prediction.
    x = x + K * y;      

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    % Update state covariance with new knowledge.
    I = eye(length(P));
    P = (I - K * H) * P;  

    % Populate the orientation quaternion matrix.
    q(i, :) = [x(1), x(2), x(3), x(4)];
    % Populate the gyroscope bias matrix.
    wb(i, :) = [x(5), x(6), x(7)];
    
end

end
% END OF FUSIONEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [marker_mbgt, T_mbgt, marker_mbcusum, T_mbcusum, marker_fsd, ...
    T_fsd, T_fsd_expanded, marker_ltsd, T_ltsd, T_ltsd_expanded, ...
    marker_amd, T_amd, marker_amvd, T_amvd, marker_ared, T_ared, ...
    marker_shod, T_shod] = build_int_markers(input_signal, ax, ay, az, ...
    gx, gy, gz, threshold_mbgt, threshold_mbcusum, threshold_fsd, ...
    threshold_ltsd, lwin_mbgt, lwin_mbcusum, lwin_fsd, lwin_ltsd, ...
    shift_fsd, shift_ltsd, lambda, lwin_amd, threshold_amd, lwin_amvd, ...
    threshold_amvd, lwin_ared, threshold_ared, lwin_shod, threshold_shod)

% FUNCTION BUILD_INT_MARKERS computes a series of motion intensity markers
% applying the MBGT, MBCUSUM, FSD, LTSD, AMD, AMVD, ARED and SHOD
% algorithms. For more information about these algorithms check the
% following publication: 'Detection of (In)activity Periods in Human Body 
% Motion Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%
% - INPUT PARAMETERS:
%   |_ 'input_signal': Signal acting as the input for the MBGT, MBCUSUM,
%                      FSD and LTSD algorithms.
%   |_ 'ax': Calibrated acceleration measured along the X-axis.
%   |_ 'ay': Calibrated acceleration measured along the Y-axis.
%   |_ 'az': Calibrated acceleration measured along the Z-axis.
%   |_ 'gx': Calibrated angular rate measured along the X-axis.
%   |_ 'gy': Calibrated angular rate measured along the Y-axis.
%   |_ 'gz': Calibrated angular rate measured along the Z-axis.
%   |_ 'threshold_mbgt': Decision threshold of the MBGT algorithm.
%   |_ 'threshold_mbcusum': Decision threshold of the MBCUSUM algorithm.
%   |_ 'threshold_fsd': Decision threshold of the FSD algorithm.
%   |_ 'threshold_ltsd': Decision threshold of the LTSD algorithm.
%   |_ 'lwin_mbgt': Size of the sliding window of the MBGT algorithm.
%   |_ 'lwin_mbcusum': Size of the sliding window of the MBCUSUM algorithm.
%   |_ 'lwin_fsd': Size of the sliding window of the FSD algorithm.
%   |_ 'lwin_ltsd': Size of the sliding window of the LTSD algorithm.
%   |_ 'shif_fsd': Overlapping of the sliding window of the FSD algorithm.
%   |_ 'shif_ltsd': Overlapping of the sliding window of the LTSD 
%                   algorithm.
%   |_ 'lambda': Normalization factor of the MBCUSUM algorithm.
%   |_ 'lwin_amd': Size of the sliding window of the AMD algorithm.
%   |_ 'threshold_amd': Decision threshold of the AMD algorithm.
%   |_ 'lwin_amvd': Size of the sliding window of the AMVD algorithm.
%   |_ 'threshold_amvd': Decision threshold of the AMVD algorithm.
%   |_ 'lwin_ared': Size of the sliding window of the ARED algorithm.
%   |_ 'threshold_ared': Decision threshold of the ARED algorithm.
%   |_ 'lwin_shod': Size of the sliding window of the SHOD algorithm.
%   |_ 'threshold_shod': Decision threshold of the SHOD algorithm.

% - OUTPUT PARAMETERS:
%   |_ 'marker_mbgt': Motion intensity marker of the MBGT algorithm.
%   |_ 'T_mbgt': Decision signal of the MBGT algorithm to which the
%                threshold is compared.
%   |_ 'marker_mbcusum': Motion intensity marker of the MBCUSUM algorithm.
%   |_ 'T_mbcusum': Decision signal of the MBCUSUM algorithm to which the
%                   threshold is compared.
%   |_ 'marker_fsd': Motion intensity marker of the FSD algorithm.
%   |_ 'T_fsd': Decision signal of the FSD algorithm to which the threshold
%               is compared.
%   |_ 'T_fsd_expanded': Decision signal of the FSD algorithm to which the
%                        threshold is compared. This signal has the same 
%                        length of the input signal. 
%   |_ 'marker_ltsd': Motion intensity marker of the LTSD algorithm.
%   |_ 'T_ltsd': Decision signal of the LTSD algorithm to which the 
%                threshold is compared.
%   |_ 'T_ltsd_expanded': Decision signal of the LTSD algorithm to which 
%                         the threshold is compared. This signal has the 
%                         same length of the input signal. 
%   |_ 'marker_amd': Motion intensity marker of the AMD algorithm.
%   |_ 'T_amd': Decision signal of the AMD algorithm to which the
%                threshold is compared.
%   |_ 'marker_amvd': Motion intensity marker of the AMVD algorithm.
%   |_ 'T_amvd': Decision signal of the AMVD algorithm to which the
%                threshold is compared.
%   |_ 'marker_ared': Motion intensity marker of the ARED algorithm.
%   |_ 'T_ared': Decision signal of the ARED algorithm to which the
%                threshold is compared.
%   |_ 'marker_shod': Motion intensity marker of the SHOD algorithm.
%   |_ 'T_shod': Decision signal of the SHOD algorithm to which the
%                threshold is compared.
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
%  Copyright (c) 2015, Alberto Olivares. All rights reserved. Based on the 
%  code in: http://goo.gl/Yww9Bb
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

% 1) Get the decision signal of the MBGT algorithm and the marker.
T_mbgt = mbgt(input_signal, lwin_mbgt);
marker_mbgt = zeros(1, length(T_mbgt));
for i = 1 : length(T_mbgt)
    if T_mbgt(i) > threshold_mbgt;
       marker_mbgt(i) = 1;
    end
end

% 2) Get the decision signal of the MBCUSUM algorithm and the marker.
T_mbcusum = mbcusum(input_signal, lwin_mbcusum, lambda);
marker_mbcusum = zeros(1, length(T_mbcusum));
for i = 1 : length(T_mbcusum)
    if T_mbcusum(i) > threshold_mbcusum;
       marker_mbcusum(i) = 1;
    end
end

% 3) Get the decision signal of the FSD algorithm and the marker.
[V_fsd, T_fsd] = fsd(input_signal, lwin_fsd, shift_fsd, 512, ...
    threshold_fsd);
[marker_fsd, T_fsd_expanded] = compEstMark(V_fsd, T_fsd, input_signal, ...
    lwin_fsd, shift_fsd);

% 4) Get the decision signal of the LTSD algorithm and the marker.
[V_ltsd, T_ltsd] = ltsd(input_signal, lwin_ltsd, shift_ltsd, 512, ...
    threshold_ltsd);
[marker_ltsd, T_ltsd_expanded] = compEstMark(V_ltsd, T_ltsd, ...
    input_signal, lwin_ltsd, shift_ltsd);

% 5) Get the decision signal of the ARED algorithm and the marker.
T_ared = ared(gx, gy, gz, lwin_ared);
marker_ared = zeros(1, length(T_ared));
for i = 1 : length(T_ared)
    if T_ared(i) > threshold_ared;
       marker_ared(i) = 1;
    end
end

% 6) Get the decision signal of the AMD algorithm and the marker.
T_amd = amd2(ax, ay, az, lwin_amd);
marker_amd = zeros(1, length(T_amd));
for i = 1 : length(T_amd)
    if T_ared(i) > threshold_amd;
       marker_amd(i) = 1;
    end
end

% 7) Get the decision signal of the AMVD algorithm and the marker.
T_amvd = amvd(ax, ay, az, lwin_amvd);
marker_amvd = zeros(1, length(T_amvd));
for i = 1 : length(T_amvd)
    if T_ared(i) > threshold_amvd;
       marker_amvd(i) = 1;
    end
end

% 8) Get the decision signal of the SHOD algorithm and the marker.
T_shod = shod(ax, ay, az, gx, gy, gz, lwin_shod);
marker_shod = zeros(1, length(T_shod));
for i = 1 : length(T_shod)
    if T_ared(i) > threshold_shod;
       marker_shod(i) = 1;
    end
end

end
% END OF BUILD_INT_MARKERS FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = mbgt(signal, lwin)

% FUNCTION MBGT applies the Memory Based Graph Theoretic Detector to 
% inertial signals in order to get a profile of the intensity of the
% motion.
%
% - INPUT PARAMETERS: 
%   |_ 'signal': Signal containing the inertial information. This signal
%                could be the magnitude of the acceleration, angular 
%                velocity, magnetic field or a linear combination of them.
%   |_ 'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_ 'T': Decision function.
%
% For more information about the MBGT algorithm see: Daniel Nikovski and 
% Ankur Jain, "Memory-Based Algorithms for Abrupt Change Detection in 
% Sensor Data Streams",Industrial Informatics, 2007 5th IEEE International 
% Conference on , vol.1, no., pp.547-552, 23-27 June 2007 .
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

L = length(signal);
for z = 1 : L - lwin 
    frame = signal(z : z + lwin);
    N = length(frame);
    Cprev = zeros(1, N + 1);
    maxim = 0;
    for i = 1 : N - 1
        beta = zeros(1, N + 1);
        for j = i + 1 : N + i - 1
            if N - i < N - j + i + 1
                beta(N - j + i + 1) = beta(N - j + i + 2) + ...
                    sqrt((frame(N - i) - frame(N - j + i + 1)) ^ 2);
            end
        end
        C = Cprev + beta;
        Cprev = C;
        frame_max = max(C);
        if frame_max > maxim
            maxim = frame_max;
        end
    end
    T(z + lwin) = maxim;
end

end
% END OF MBGT FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = mbcusum(signal, lwin, lambda)

% FUNCTION MBCUSUM applies the Memory Based CUSUM Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'signal': Signal containing the inertial information. This signal
%                could be the magnitude of the acceleration, angular 
%                velocity, magnetic field or a linear combination of them.
%   |_ 'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_ 'T': Decision function.
%
% For more information about the FSD algorithm see.
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

L = length(signal);
for z = 1 : L - lwin 
    frame = signal(z : z + lwin - 1);
    N = length(frame);
    vprev = zeros(1, N + 1);
    mu = zeros(1, N + 1);
    v = zeros(1, N + 1);
    s = zeros(1, N + 1);
    maxim = 0;
    for i = N - 1 : -1 : 1
        for j = N : -1 : i + 1
            sum_weights = 0;
            if i < j
                for k = j : N
                    value = 1 / (N * (2 * lambda ^ 2 * pi) ^ (1 /2 )) ...
                        * exp(-1 / 2 * (abs(frame(j) - frame(k)) / ...
                        lambda) ^ 2);
                    sum_weights = sum_weights + value;
                end
            end
            mu(j) = sum_weights;
            v(j) = vprev(j) + 1 / (N * (2 * lambda ^ 2 * pi) ^ (1 / 2)) ...
                * exp(-1 / 2 * (abs(frame(j) - frame(i)) / lambda) ^ 2);
            vprev = v;
            s(j) = s(j + 1) + log(mu(j)) - log(v(j)) + log(j - i) - ...
                log(N - j + 1);
        end
        frame_max = max(s);
        if frame_max > maxim
            maxim = frame_max;
        end
    end
     T(z) = maxim;
end

end
% END OF MBCUSUM FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [V, T] = fsd(x, LWIN, SHIFT, NFFT, Threshold)

% FUNCTION FSD applies the Framed Spectrum Detector to inertial signals in 
% order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'x': Signal containing the inertial information. This signal could 
%           be the magnitude of the acceleration, angular velocity, 
%           magnetic field or a linear combination of them.
%   |_ 'LWIN': Length of sliding window.
%   |_ 'SHIFT': Overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%   |_ 'Threshold': Decision threshold.
% 
% - OUTPUT PARAMETERS:
%   |_ 'V': Reduced activity marker.
%   |_ 'T': Decision function.
%
% For more information about the FSD algorithm check the following 
% publication: 'Detection of (In)activity Periods in Human Body Motion 
% Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Javier Ramírez Pérez de Inestrosa (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  javierrp@ugr.es
%   |_ Dr. Dr. Juan Manuel Górriz Sáez (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  gorriz@ugr.es
%   |_ Dr. Alberto Olivares (Reorganization):
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
%  Copyright (c) 2015, Javier Ramírez, Juan Manuel Górriz and 
%  Alberto Olivares. All rights reserved.
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

% Computation of the number of frames.
Nframes = floor((length(x) - LWIN) / SHIFT) + 1;

% Computation of spectrum.
SP = spectrum2(x, LWIN, SHIFT, NFFT);

% Initialization of the reduced marker.
V = zeros(Nframes, 1);

% Noise averaging initizalization time (frames)
FI = 1;

% Noise estimation.
NE = 0;
alfa = 0.98;

% Initialization of the decision function.
T = zeros(Nframes, 1);

for frame = 1 : Nframes
       
    if (frame <= FI)
       NE = (1 - 1 /frame) * NE + 1 / frame * SP(:, frame);
       V(frame) = 0;
    else
       T(frame) = 10 * log10(1 / (NFFT / 2) * sum(SP(:, frame) .^ 2 ./ ...
           NE .^ 2));
       if (T(frame) > Threshold)
           V(frame) = 1; 
       else
           V(frame) = 0;
           NE = alfa * NE + (1 - alfa) * SP(:, frame);
       end
    end
end

end
% END OF FSD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [y1, y2] = compEstMark(V, T, x, LWIN, SHIFT)

% FUNCTION COMPESTMARK transforms the FSD and LTSD decision signals into a
% marker of the same length of the input signal.
%
% - INPUT PARAMETERS:
%   |_'V': Decision function computed with LTSD or FSD.
%   |_'T': Marker of reduced length.
%   |_'x': Signal used as input to compute the figure of merit with LTSD or
%          FSD, i.e. acceleration or angular rate magnitude, sum of 
%          magnitudes or product of magnitudes.
%   |_'LWIN': Window length used to compute the figure of merit with LTSD 
%             or FSD.
%   |_'SHIFT': Overlapping used to compute the figure of merit with LTSD or
%              FSD.
% 
% - OUTPUT PARAMETERS:
%   |_'y1': Intensity binary marker (extended length).
%   |_'y2': Decision function (extended length).
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Javier Ramírez Pérez de Inestrosa (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  javierrp@ugr.es
%   |_ Dr. Dr. Juan Manuel Górriz Sáez (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  gorriz@ugr.es
%   |_ Dr. Alberto Olivares (Reorganization):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/22/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Javier Ramírez, Juan Manuel Górriz and 
%  Alberto Olivares. All rights reserved.
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

y1 = zeros(1, length(x));
y2 = zeros(1, length(x));

N = min(floor((length(x) - LWIN) / SHIFT) + 1, length(V));

for i = 1 : N
    y2((i - 1) * SHIFT + (1 : LWIN)) = T(i);
   if (V(i) == 0)
          y1((i - 1) * SHIFT + (1 : LWIN)) = 0; 
   else
          y1((i - 1) * SHIFT + (1 : LWIN)) = 1; 
   end
end

end
% END OF COMPESTMARK FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [V, T] = ltsd(x, LWIN, SHIFT, NFFT, threshold)

% FUNCTION LTSD applies the Long Term Spectrum Detector to inertial signals
% in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS: 
%   |_ 'x': Signal containing the inertial information. This signal could 
%           be the magnitude of the acceleration, angular velocity, 
%           magnetic field or a linear combination of them.
%   |_ 'LWIN': Length of sliding window.
%   |_ 'SHIFT': Overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%   |_ 'threshold': Decision threshold.
% 
% - OUTPUT PARAMETERS:
%   |_ 'V': Reduced activity marker.
%   |_ 'T': Decision function.
%
% For more information about the LTSD algorithm check the following 
% publication: 'Detection of (In)activity Periods in Human Body Motion 
% Using Inertial Sensors: A Comparative Study' by A. Olivares et al.
% The article can be accessed in the following URL: 
% http://www.mdpi.com/1424-8220/12/5/5791
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Javier Ramírez Pérez de Inestrosa (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  javierrp@ugr.es
%   |_ Dr. Dr. Juan Manuel Górriz Sáez (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  gorriz@ugr.es
%   |_ Dr. Alberto Olivares (Reorganization):
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
%  Copyright (c) 2015, Javier Ramírez, Juan Manuel Górriz and 
%  Alberto Olivares. All rights reserved.
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

% No-delay backward long-term spectral envelope.
N = 8;

Nframes = floor((length(x) - LWIN) / SHIFT) + 1;
SP = spectrum2(x, LWIN, SHIFT, NFFT);
V = zeros(Nframes, 1);

% Noise averaging initizalization time (frames).
FI = 1;

% Noise estimation.
NE = 0;
alfa = 0.98;

% Decision function.
T = zeros(Nframes, 1);
for frame = 1 : Nframes
       
    if (frame <= FI)
       NE = (1 - 1 / frame) * NE + 1 / frame * SP(:, frame);
       V(frame) = 0;
    else
       N1 = min(frame - 1, N);
       LTSE = max(SP(:, frame - N1 : frame)')'; 
          
       T(frame) = 10 * log10(1 / (NFFT / 2) * sum(LTSE .^ 2 ./ NE .^ 2));
       if (T(frame) > threshold)
           V(frame) = 1; 
       else
           V(frame) = 0;
           NE = alfa * NE + (1 - alfa) * SP(:, frame);
       end
    end
end

end
% END OF LTSD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = amd2(ax, ay, az, lwin)

% FUNCTION AMD2 applies the Acceleration Magnitude Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the AMD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

% Set noise variance of accelerometer and gyroscope.
var_a = 0.009;

% Compute AMD algorithm.
for i = 1 : length(ax) - lwin     
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);
    N = length(ax_frame);
    mod_aF = sqrt(ax_frame .^ 2 + ay_frame .^ 2 + az_frame .^ 2);
    T(i : i + lwin) = 1 / (N * var_a) * (mod_aF - 1) .^ 2;  
end

end
% END OF AMD2 FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = amvd(ax, ay, az, lwin)

% FUNCTION AMVD applies the Acceleration Magnitude Variance Detector to 
% inertial signals in order to get a profile of the intensity of the 
% motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the AMVD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

% Set variance of accelerometer noise.
var_a = 0.009;

% Compute AMVD algorithm.
for i = 1 : length(ax) - lwin   
    
    % Build the signal window.
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);

    % Compute the average acceleration values of the window.
    av_x = mean(ax_frame);
    av_y = mean(ay_frame);
    av_z = mean(az_frame);
    N = length(ax_frame);

    % Subtract the average acceleration to each sample of the window.
    axC = ax_frame - av_x;
    ayC = ay_frame - av_y;
    azC = az_frame - av_z;
    
    % Compute the magnitude of the modified acceleration.
    mod_aC = sqrt(axC .^ 2 + ayC .^ 2 + azC .^ 2);

    % Compute the detection signal.
    T(i : i + lwin) = 1 / (N * var_a) * mod_aC .^ 2;   
end

end
% END OF AMVD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = ared(gx, gy, gz, lwin)

% FUNCTION ARED applies the Angular Rate Energy Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'gx': Angular rate signal (X axis).
%   |_'gy': Angular rate signal (Y axis).
%   |_'gz': Angular rate signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the ARED algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

% Compute ARED algorithm.
for i = 1 : length(gx) - lwin   
    gx_frame = gx(i : i + lwin);
    gy_frame = gy(i : i + lwin);
    gz_frame = gz(i : i + lwin);
    N = length(gx_frame);
    mod_gF = sqrt(gx_frame .^ 2 + gy_frame .^ 2 + gz_frame .^ 2);
    T(i : i + lwin) = 1 / N *(mod_gF) .^ 2;  
end

end
% END OF ARED FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function T = shod(ax, ay, az, gx, gy, gz, lwin)

% FUNCTION SHOD applies the Stance Hypothesis Optimal Detector to inertial 
% signals in order to get a profile of the intensity of the motion.
%
% - INPUT PARAMETERS:
%   |_'ax': Acceleration signal (X axis).
%   |_'ay': Acceleration signal (Y axis).
%   |_'az': Acceleration signal (Z axis).
%   |_'gx': Angular rate signal (X axis).
%   |_'gy': Angular rate signal (Y axis).
%   |_'gz': Angular rate signal (Z axis).
%   |_'lwin': Length of sliding window.
% 
% - OUTPUT PARAMETERS:
%   |_'T': Decision function.
%
% For more information about the SHOD algorithm see: I Skog et al. 
% "Zero-Velocity Detection-An algorithm Evaluation", Biomedical 
% Engineering, IEEE Transactions obn., vol.57, no.11, 2010.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

% Set noise variance of accelerometer and gyroscope.
var_a = 0.009;
var_g = 5;

% Compute SHOD algorithm.
for i = 1 : length(ax) - lwin   
    
    gx_frame = gx(i : i + lwin);
    gy_frame = gy(i : i + lwin);
    gz_frame = gz(i : i + lwin);
    
    ax_frame = ax(i : i + lwin);
    ay_frame = ay(i : i + lwin);
    az_frame = az(i : i + lwin);
    N = length(ax_frame);
    
    av_x = mean(ax_frame);
    av_y = mean(ay_frame);
    av_z = mean(az_frame);
    mod_av_a = sqrt(av_x ^ 2 + av_y ^ 2 + av_z ^ 2);
    term_ax = ax_frame - 1 * (av_x / mod_av_a);
    term_ay = ay_frame - 1 * (av_y / mod_av_a);
    term_az = az_frame - 1 * (av_z / mod_av_a);
    mod_a = sqrt(term_ax .^ 2 + term_ay .^ 2 + term_az .^ 2) .^ 2;
    
    mod_g = sqrt(gx_frame .^ 2 + gy_frame .^ 2 + gz_frame .^ 2) .^ 2;
    
    T(i : i + lwin) = 1 / N * (1 / var_a * mod_a + 1 / var_g * mod_g);  
    
end

end
% END OF SHOD FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function SP = spectrum2(x, LWIN, SHIFT, NFFT)

% FUNCTION SPECTRUM computes the spectrum of a given signal. 
%
% - INPUT PARAMETERS:
%   |_ 'x': Input signal. 
%   |_ 'LWIN': Length of the sliding window.
%   |_ 'SHIFT': Length of the overlapping of the sliding window.
%   |_ 'NFFT': Resolution of the FFT.
%
% - OUTPUT PARAMETERS:
%   |_ 'SP': Spectrum of the signal.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Dr. Javier Ramírez Pérez de Inestrosa (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  javierrp@ugr.es
%   |_ Dr. Dr. Juan Manuel Górriz Sáez (Original version):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  gorriz@ugr.es
%   |_ Dr. Alberto Olivares (Reorganization):
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
%  Copyright (c) 2015, Javier Ramírez, Juan Manuel Górriz and 
%  Alberto Olivares. All rights reserved.
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

Nframes = floor((length(x) - LWIN) / SHIFT) + 1;

SP = zeros(NFFT / 2, Nframes);

for frame = 1 : Nframes
    r = (frame - 1) * SHIFT + [1 : LWIN];
    s_w = hamming(LWIN) .* x(r);   
    magnitude = abs(fft(s_w, NFFT));
    SP(:, frame) = magnitude(1 : NFFT / 2);
end

end
% END OF SPECTRUM FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeGKF(observation, ang_vel, frequency,...
    ref_angle, p0, rmse_off, chosen_marker)

% FUNCTION OPTIMIZEGKF uses the ANMS algorithm to find the optimal 
% parameters of the Gated Kalman Filter which minimize the error between
% the actual orientation angle and the one estimated with GKF. 
% 
% - INPUT PARAMETERS:
%   |_ 'observation': Observation of the Kalman Filter.
%   |_ 'ang_vel': Angular velocity (measured with the gyroscope).
%   |_ 'frequency': Sampling frequency.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'rmse_off': Offset in the RMSE computation.
%   |_ 'chosen_marker': Intensity marker.
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

% 1) Set variables.
% -------------------------------------------------------------------------
global gyro;        gyro = ang_vel;
global obs;         obs = observation;
global true_angle;  true_angle = ref_angle;
global freq;        freq = frequency;
global rmse_offset; rmse_offset = rmse_off;
global marker;      marker = chosen_marker;      

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp(['Optimizing parameters of Gated Kalman Filter',...
    ' (it may take a while) ...']);

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin, fmin, ct] = ANMS(@eofGKF, p0, tol, max_feval);

end
% END OF OPTIMIZEGKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofGKF(p)

% FUNCTION EOFKGKF Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the Gated Kalman Filter. 
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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
global gyro;
global obs;
global freq;
global true_angle;
global rmse_offset;
global marker;

% Extract the first parameter to be minimized (measurement noise variance 
% gain). 
alpha1 = p(1);
alpha2 = p(2);
beta1 = p(3);
beta2 = p(4);

% Compute the orientation angle estimation. 
angle_GKF = fusionGKF(gyro, obs, freq, alpha1, alpha2, beta1, beta2, ...
    marker);

% Define error function.
F = compute_rmse(true_angle, 180 / pi * angle_GKF, rmse_offset);

end
% END OF EOFGKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [x1, x2, v_output, K_output] = fusionGKF(gyro, obs, freq, ...
    alpha1, alpha2, beta1, beta2, marker)

% FUNCTION FUSIONGKF Applies a Gated Kalman Filter sensor fusion approach 
% to estimate the orientation of a body using the acceleration, magnetic 
% field and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'gyro': vector containing the gyroscope signal for a determined 
%               axis. 
%   |_ 'obs': vector containing the angle values obtained via the linear
%             acceleration relations.
%   |_ 'freq': sampling frecuency. Must be real positive.
%   |_ 'alpha1': Weighing factor which multiplies the variance of the
%                measurement noise when the intensity of the motion is low. 
%                It is a parameter which tunes the filter.
%   |_ 'alpha2': Weighing factor which multiplies the variance of the
%                measurement noise when the intensity of the motion is 
%                high. It is a parameter which tunes the filter.
%   |_ 'beta1': Weighing factor which multiplies the variance of the 
%               process noise when the intensity of the motion is low. It 
%               is a parameter which tunes the filter.
%   |_ 'beta2': Weighing factor which multiplies the variance of the 
%               process noise when the intensity of the motion is high. It 
%               is a parameter which tunes the filter.  
%
% - OUTPUT PARAMETERS:
%   |_ 'x1': estimated orientation angle. (First element of the state
%            vector).
%   |_ 'x2': estimated gyroscope bias. (Second element of the state
%            vector).
%   |_ 'v_output': Vector containing the state vector corrected 'a
%                  posteriori' by the observation.
%   |_ 'K_output': Vector containing the value of the Kalman gain for each
%                  time instant.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

% 1) Definition of variables.
% -------------------------------------------------------------------------
obsVar = var(obs);
gyroVar = var(gyro);
len = length(gyro);

% Sampling period.
dt = 1 / freq;

% Initialization of covariance matrix.
P = [1 0; 0 1];

% Set measurement variance.
Rk1 = alpha1 * obsVar;
Rk2 = alpha2 * obsVar;

% Assuming process noise is white for each component, covariance matrix of
% process noise must be a diagonal matrix, with noise variance for each
% component in each position of the diagonal.
Q1 = beta1 * [gyroVar 0; 0 obsVar];
Q2 = beta2 * [gyroVar 0; 0 obsVar];

% Set state transition matrix.
A = [1 -dt; 0 1];

% Set measurement matrix. 
C = [1 0];

% Initialize state vector.
X = [0; 0];
x1 = zeros(len, 1);
x2 = zeros(len, 1);

% 2) Body of the Kalman Filter.
% -------------------------------------------------------------------------
for i = 1 : len
    
    if marker(i) == 0
        Rk = Rk1;
        Q = Q1;
    end
    if marker(i) == 1
        Rk = Rk2;
        Q = Q2;
    end
    
    % Prediction phase.
    X(1) = X(1) + (gyro(i) - X(2))*dt;
    X(2) = X(2);
    P = A * P * A' + Q;

    % Update phase.
    v = obs(i) - C * X;
    v_output(i) = v;
    Sk = C * P * C' + Rk;
    K = (P * C') / Sk;
    K_output(i, :) = K';
    X = X + K * v;
    P = P - K * Sk * K';
    x1(i) = X(1);
    x2(i) = X(2);   
end

end
% END OF FUSIONGKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeGEKF(axC, ayC, azC, gxC, gyC, gzC, ...
    frequency, ref_angle, p0, rmse_off, selected_angle, chosen_marker)

% FUNCTION OPTIMIZEEKF uses the ANMS algorithm to find the optimal 
% parameters of the Gated Extended Kalman Filter which minimize the error 
% between the actual orientation angle and the one estimated with GEKF. 
% 
% - INPUT PARAMETERS:
%   |_ 'axC': Calibrated acceleration measured along the X-axis.
%   |_ 'ayC': Calibrated acceleration measured along the Y-axis.
%   |_ 'azC': Calibrated acceleration measured along the Z-axis.
%   |_ 'gxC': Calibrated angular rate measured along the X-axis.
%   |_ 'gyC': Calibrated angular rate measured along the Y-axis.
%   |_ 'gzC': Calibrated angular rate measured along the Z-axis.
%   |_ 'frequency': Sampling frequency.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'rmse_off': Offset in the RMSE computation.
%   |_ 'selected_angle': Orientation angle being estimated. Possible values
%                        are 'pitch', 'roll' and 'yaw'.
%   |_ 'chosen_marker': Intensity marker.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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
global gyro;        gyro = [gxC gyC gzC];
global acc;         acc = [azC ayC axC];
global freq;        freq = frequency;
global true_angle;  true_angle = ref_angle;
global rmse_offset; rmse_offset = rmse_off;
global chosen_angle;chosen_angle = selected_angle;
global marker;      marker = chosen_marker;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp(['Optimizing parameters of Gated Extended Kalman Filter',...
    ' (it may take a while)...']);

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-5;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin,fmin,ct] = ANMS(@eofGEKF, p0, tol, max_feval);
 
end
% END OF OPTIMIZEGEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofGEKF(p)

% FUNCTION EOFGEKF Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the Gated Extended Kalman Filter. 
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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
global gyro;
global acc;
global true_angle;
global freq;
global chosen_angle;
global marker;
global rmse_offset;

% Extract the value of the first parameter (variance of gyroscope bias).
sigma_w1 = p(1);
sigma_w2 = p(3);

% Extract the value of the second parameter (variance of the observation).
sigma_obs1 = p(2);
sigma_obs2 = p(4);

% Set sampling period. 
dt = 1 / freq;

% Initialize orientation quaternion. 
ini = [1 0 0 0];

% Call the EKF routine to compute the orientation quaternion.
q = fusionGEKF(-acc, gyro, dt, sigma_w1, sigma_w2, sigma_obs1, ...
    sigma_obs2, ini, marker);

% Transform the orientation quaternion to Euler angles.
roll_GEKF = zeros(1, length(q));
pitch_GEKF = zeros(1, length(q));
yaw_GEKF = zeros(1, length(q));
for i = 1 : length(q)
     eulerAngles = quat_to_euler(q(i, :));
     roll_GEKF(i) = eulerAngles(1);
     pitch_GEKF(i) = eulerAngles(2);
     yaw_GEKF(i) = eulerAngles(3);
end

% Select the chosen angle. 
if strcmp(chosen_angle, 'pitch')
    angle_GEKF = pitch_GEKF;
elseif strcmp(chosen_angle, 'roll')
    angle_GEKF = roll_GEKF;
elseif strcmp(chosen_angle, 'yaw')
    angle_GEKF = yaw_GEKF;
end

% Define the error function to be minimized. 
F = compute_rmse(true_angle, 180 / pi * angle_GEKF', rmse_offset);
end
% END OF EOFGEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [q, p1, p2 , p3, p4, p5, p6, p7, wb] = fusionGEKF(a, w, dt, ...
    sigma_w1, sigma_w2, sigma_obs1, sigma_obs2, quat_ini, marker)

% FUNCTION FUSIONGEKF Applies a Kalman Filter sensor fusion approach to
% estimate the orientation of a body using the acceleration, magnetic field
% and angular rate measured with a triaxial accelerometer, a triaxial
% magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'a': Acceleration matrix.
%   |_ 'w': Angular velocity matrix.
%   |_ 'dt': Sampling period.
%   |_ 'sigma_w1': Variance of the process noise. This is a tuning 
%                  parameter of the filter when the motion intensity is
%                  low.
%   |_ 'sigma_w2': Variance of the process noise. This is a tuning 
%                  parameter of the filter when the motion intensity is
%                  high.
%   |_ 'sigma_obs1': Variance of the measurement noise. This is a tuning
%                    parameter of the filter when the motion intensity is
%                    low.
%   |_ 'sigma_obs2': Variance of the measurement noise. This is a tuning
%                    parameter of the filter when the motion intensity is
%                    high.
%   |_ 'quat_ini': Quaternion determining the initial orientation.
%   |_ 'marker': Intensity marker.
%
% - OUTPUT PARAMETERS:
%   |_ 'q': Estimated orientation quaternion.
%   |_ 'p1': Element (1,1) of the covariance matrix.
%   |_ 'p2': Element (2,2) of the covariance matrix.
%   |_ 'p3': Element (3,3) of the covariance matrix.
%   |_ 'p4': Element (4,4) of the covariance matrix.
%   |_ 'p5': Element (5,5) of the covariance matrix.
%   |_ 'p6': Element (6,6) of the covariance matrix.
%   |_ 'p7': Element (7,7) of the covariance matrix.
%   |_ 'wb': Vector containing the bias of the angular velocities.
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
%  Copyright (c) 2015, Alberto Olivares. All rights reserved. Based on the 
%  code in: http://goo.gl/Yww9Bb
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
% Initialize orientation quaternion.
q = zeros(length(a(:, 1)), 4);

% Initialize gyroscope bias vector. 
wb = zeros(length(a(:, 1)), 3);

% Initialize the system covariance matrix.
P = eye(7) * 10000; 

% Initialize the state vector
x = [quat_ini 0 0 0]';

% 2) Body of the EKF.
% -------------------------------------------------------------------------
for i = 1 : length(a(:, 1))
    
    if marker(i) == 0
        % Initialize the process noise covariance matrix.
        Q = [0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, sigma_w1,   0,          0;
             0, 0, 0, 0, 0,         sigma_w1,    0;
             0, 0, 0, 0, 0,         0,          sigma_w1];

        % Initialize the measurement noise covariance matrix.
        R = [sigma_obs1  0           0  ;
             0          sigma_obs1   0  ;
             0          0           sigma_obs1;];

    elseif marker(i) == 1
            % Initialize the process noise covariance matrix.
        Q = [0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, 0,         0,          0;
             0, 0, 0, 0, sigma_w2,   0,          0;
             0, 0, 0, 0, 0,         sigma_w2,    0;
             0, 0, 0, 0, 0,         0,          sigma_w2];

        % Initialize the measurement noise covariance matrix.
        R = [sigma_obs2  0           0  ;
             0          sigma_obs2   0  ;
             0          0           sigma_obs2;];

    end
   
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    wxb = x(5);
    wyb = x(6);
    wzb = x(7);
    wx = w(i, 1);
    wy = w(i, 2);
    wz = w(i, 3);
    clear z
    z(1) = a(i, 1);
    z(2) = a(i, 2);
    z(3) = a(i, 3);
    z = z';

    % Populate F jacobian
    F = [              1,  (dt/2)*(wx-wxb),  (dt/2)*(wy-wyb), ...
         (dt/2)*(wz-wzb),       -(dt/2)*q1,       -(dt/2)*q2, -(dt/2)*q3;
        -(dt/2)*(wx-wxb),                1,  (dt/2)*(wz-wzb), ...
        -(dt/2)*(wy-wyb),        (dt/2)*q0,        (dt/2)*q3, -(dt/2)*q2;
        -(dt/2)*(wy-wyb), -(dt/2)*(wz-wzb),                1, ...
         (dt/2)*(wx-wxb),       -(dt/2)*q3,        (dt/2)*q0, (dt/2)*q1;
        -(dt/2)*(wz-wzb),  (dt/2)*(wy-wyb), -(dt/2)*(wx-wxb), ...  
                       1,        (dt/2)*q2,       -(dt/2)*q1, (dt/2)*q0;
                       0,                0,                0, ...
                       0,                1,                0,          0;
                       0,                0,                0, ...
                       0,                0,                1,          0;
                       0,                0,                0, ...
                       0,                0,                0,          1;];
                   
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Predicted state estimate, x = f(x,u).
    x = [q0 - (dt/2) * (-q1*(wx-wxb) - q2*(wy-wyb) - q3*(wz-wzb));
         q1 - (dt/2) * ( q0*(wx-wxb) + q3*(wy-wyb) - q2*(wz-wzb));
         q2 - (dt/2) * (-q3*(wx-wxb) + q0*(wy-wyb) + q1*(wz-wzb));
         q3 - (dt/2) * ( q2*(wx-wxb) - q1*(wy-wyb) + q0*(wz-wzb));
         wxb;
         wyb;
         wzb;];

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);

    % Predicted covariance estimate.
    P = F * P * F' + Q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Normalize accelerometer and magnetometer measurements.
    acc_norm = sqrt(z(1) ^ 2 + z(2) ^ 2 + z(3) ^ 2);
    z(1) = z(1) / acc_norm;
    z(2) = z(2) / acc_norm;
    z(3) = z(3) / acc_norm;

    h = [-2 * (q1 * q3 - q0 * q2);
         -2 * (q2 * q3 + q0 * q1);
         -2 * q0 ^ 2 + 1 - 2 * q3 ^ 2];

    % Measurement residual: y = z - h(x), where h(x) is the matrix that 
    % maps the state onto the measurement.
    y = z - h;

    % The H matrix maps the measurement to the states.
    H = [ 2 * q2, -2 * q3,  2 * q0, -2 * q1, 0, 0, 0;
         -2 * q1, -2 * q0, -2 * q3, -2 * q2, 0, 0, 0;
         -4 * q0,       0,       0, -4 * q3, 0, 0, 0;]; 
     
    % Measurement covariance update.
    S = H * P * H' + R;

    % Calculate Kalman gain.
    K = P * H'/ S;

    % Corrected model prediction.
    x = x + K * y;      

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    % Update state covariance with new knowledge.
    I = eye(length(P));
    P = (I - K * H) * P;  

    % Populate the orientation quaternion matrix.
    q(i,:) = [x(1), x(2), x(3), x(4)];
    
    % Populate the gyroscope bias matrix.
    wb(i,:) = [x(5), x(6), x(7)];
    
    % Get the diagonal elements of the covariance matrix.
    p1(i) = P(1, 1);
    p2(i) = P(2, 2);
    p3(i) = P(3, 3);
    p4(i) = P(4, 4);
    p5(i) = P(5, 5);
    p6(i) = P(6, 6);
    p7(i) = P(7, 7);
end

end
% END OF FUSIONGEKF FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [xmin, fmin, ct] = optimizeGEKF_QUEST(angular_vel, obs, p0, ...
    ref_angle, selected_angle, frequency, rmse_off, chosen_marker)

% FUNCTION OPTIMIZEGEKF_QUEST uses the ANMS algorithm to find the optimal 
% parameters of the Extended Kalman Filter + QUEST which minimize the error
% between the actual orientation angle and the one estimated with EKF. 
% 
% - INPUT PARAMETERS:
%   |_ 'angular_vel': Angular velocity (measured with the gyroscope).
%   |_ 'obs': Observation of the Kalman Filter.
%   |_ 'p0': Initial value of the parameters to be optimized.
%   |_ 'ref_angle': Orientation angle reference measured with the
%                   mechanical device.
%   |_ 'selected_angle': Orientation angle being estimated. Possible values
%                        are 'pitch', 'roll' and 'yaw'.
%   |_ 'frequency': Sampling frequency.
%   |_ 'rmse_off': Offset in the RMSE computation.
%   |_ 'chosen_marker': Intensity marker.
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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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
global gyro;            gyro = angular_vel;
global quat_obs;        quat_obs = obs;
global true_angle;      true_angle = ref_angle;
global chosen_angle;    chosen_angle = selected_angle;
global freq;            freq = frequency;
global rmse_offset;     rmse_offset = rmse_off;
global marker;          marker = chosen_marker;

% 2) Call the minimization routine.
% -------------------------------------------------------------------------
disp(['Optimizing parameters of Gated Extended Kalman Filter (QUEST)',...
    ' (it may take a while) ...']);

% Set tolerance limit between the minimum values found in subsequent
% iterations of the algorithm. 
tol = 10^-6;

% Set maximum number of evaluations of the error function.
max_feval = 5000;

% Call ANMS algorithm to perform optimization (find optimal parameters of
% the QUEST algorithm which minimize the error function).
[xmin, fmin, ct] = ANMS(@eofGEKF_QUEST, p0, tol, max_feval);

end
% END OF OPTIMIZEGEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function F = eofGEKF_QUEST(p)

% FUNCTION EOFGEKF_QUEST Computes the error function to be minimized: RMSE
% between the actual angle and the estiamted orientation angle computed
% with the GEKF-QUEST Algorithm. 

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
%   |_ Date: 01/22/2015 (mm/dd/yy). 
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

global gyro;
global freq;
global quat_obs;
global true_angle;
global chosen_angle;
global marker;
global rmse_offset;

dt = 1 / freq;

sigma_w1 = p(1);
sigma_obs1 = p(2);
sigma_w2 = p(3);
sigma_obs2 = p(4);
quat_ini = quat_obs(1, :);

quat_est = fusionGEKF_QUEST(quat_obs, gyro, dt, sigma_w1, sigma_w2, ...
    sigma_obs1, sigma_obs2, quat_ini, marker);

roll_GEKF_QUEST = zeros(1, length(quat_est));
pitch_GEKF_QUEST = zeros(1, length(quat_est));
yaw_GEKF_QUEST = zeros(1, length(quat_est));
for i = 1 : length(quat_est)
     eulerAngles = quat_to_euler(quat_est(i,:));
     roll_GEKF_QUEST(i) = eulerAngles(1);
     pitch_GEKF_QUEST(i) =- eulerAngles(2);
     yaw_GEKF_QUEST(i) = eulerAngles(3);
end

% Select the chosen angle. 
if strcmp(chosen_angle, 'pitch')
    angle_GEKF_QUEST = pitch_GEKF_QUEST;
elseif strcmp(chosen_angle, 'roll')
    angle_GEKF_QUEST = roll_GEKF_QUEST;
elseif strcmp(chosen_angle, 'yaw')
    angle_GEKF_QUEST = yaw_GEKF_QUEST;
end

% Define the error function to be minimized. 
F = compute_rmse(true_angle, 180 / pi * angle_GEKF_QUEST', rmse_offset);

end
% END OF EOFGEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [q, wb, y_output] = fusionGEKF_QUEST(quat, w, dt, sigma_w1, ...
    sigma_w2, sigma_obs1, sigma_obs2, quat_ini, marker)

% FUNCTION FUSIONGEKF_QUEST Applies a Gated Kalman Filter sensor fusion 
% approach to estimate the orientation of a body using the acceleration, 
% magnetic field and angular rate measured with a triaxial accelerometer, a
% triaxial magnetometer and a triaxial gyroscope. 
%
% - INPUT PARAMETERS:
%   |_ 'quat': Quaternion containing the orientation estimated by the QUEST
%              algorithm. These estimates are used as the observation.
%   |_ 'w': Angular velocity matrix.
%   |_ 'dt': Sampling period.
%   |_ 'sigma_w1': Variance of the process noise. This is a tuning 
%                  parameter of the filter when the motion intensity is
%                  low.
%   |_ 'sigma_w2': Variance of the process noise. This is a tuning 
%                  parameter of the filter when the motion intensity is
%                  high.
%   |_ 'sigma_obs1': Variance of the measurement noise. This is a tuning
%                    parameter of the filter when the motion intensity is
%                    low.
%   |_ 'sigma_obs2': Variance of the measurement noise. This is a tuning
%                    parameter of the filter when the motion intensity is
%                    high.
%   |_ 'quat_ini': Quaternion determining the initial orientation.
%   |_ 'marker': Intensity marker.
%
% - OUTPUT PARAMETERS:
%   |_ 'q': Estimated orientation quaternion.
%   |_ 'wb': Vector containing the bias of the angular velocities.
%   |_ 'y_output': Measurement residual of every time instant.
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
%  Copyright (c) 2015, Alberto Olivares. All rights reserved. Based on the 
%  code in: http://goo.gl/Yww9Bb
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
% Initialize orientation quaternion.
q = zeros(length(quat(:, 1)), 4);

% Initialize gyroscope bias vector. 
wb=zeros(length(w(:, 1)), 3);

% Initialize the system covariance matrix.
P = eye(length(7)) * 10000; 

% Initialize the state vector
x = [quat_ini 0 0 0]';

% 2) Body of the EKF.
% -------------------------------------------------------------------------
for i = 1 : length(quat(:, 1))
    
    if marker(i) == 0
        sigma_w = sigma_w1;
        sigma_obs = sigma_obs1;
    elseif marker(i) == 1
        sigma_w = sigma_w2;
        sigma_obs = sigma_obs2;
    end
    
    % Initialize the process noise covariance matrix.
    Q = [0, 0, 0, 0, 0,         0,          0;
         0, 0, 0, 0, 0,         0,          0;
         0, 0, 0, 0, 0,         0,          0;
         0, 0, 0, 0, 0,         0,          0;
         0, 0, 0, 0, sigma_w,   0,          0;
         0, 0, 0, 0, 0,         sigma_w,    0;
         0, 0, 0, 0, 0,         0,          sigma_w];

    % Initialize the measurement noise covariance matrix.
    R = [sigma_obs  0           0          0          0        0        0;
         0          sigma_obs   0          0          0        0        0;
         0          0           sigma_obs  0          0        0        0;
         0          0           0          sigma_obs  0        0        0;
         0          0           0          0          sigma_w  0        0;
         0          0           0          0          0        sigma_w  0
         0          0           0          0          0        0        ...
                                                                  sigma_w];
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    wxb = x(5);
    wyb = x(6);
    wzb = x(7);
    wx = w(i, 1);
    wy = w(i, 2);
    wz = w(i, 3);
    
    % Populate F jacobian
    F = [              1,  (dt/2)*(wx-wxb),  (dt/2)*(wy-wyb), ...
         (dt/2)*(wz-wzb),       -(dt/2)*q1,       -(dt/2)*q2, -(dt/2)*q3;
        -(dt/2)*(wx-wxb),                1,  (dt/2)*(wz-wzb), ...
        -(dt/2)*(wy-wyb),        (dt/2)*q0,        (dt/2)*q3, -(dt/2)*q2;
        -(dt/2)*(wy-wyb), -(dt/2)*(wz-wzb),                1, ...
         (dt/2)*(wx-wxb),       -(dt/2)*q3,        (dt/2)*q0, (dt/2)*q1;
        -(dt/2)*(wz-wzb),  (dt/2)*(wy-wyb), -(dt/2)*(wx-wxb), ...  
                       1,        (dt/2)*q2,       -(dt/2)*q1, (dt/2)*q0;
                       0,                0,                0, ...
                       0,                1,                0,          0;
                       0,                0,                0, ...
                       0,                0,                1,          0;
                       0,                0,                0, ...
                       0,                0,                0,          1;];  
                   
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTION PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Predicted state estimate, x = f(x,u).
    x = [q0 - (dt / 2) * (-q1 * (wx - wxb) - q2 * (wy - wyb) - q3 * (wz ...
         - wzb));
         q1 - (dt / 2) * (q0 * (wx - wxb) + q3 * (wy - wyb) - q2 * (wz ...
         - wzb));
         q2 - (dt / 2) * (-q3 * (wx - wxb) + q0 * (wy - wyb) + q1 * (wz ...
         - wzb));
         q3 - (dt / 2) * ( q2 * (wx - wxb) - q1 * (wy - wyb) + q0 * (wz ...
         - wzb));
         wxb;
         wyb;
         wzb;];

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);

    % Predicted covariance estimate.
    P = F * P * F' + Q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h = [q0; q1; q2; q3; 0; 0; 0];

    % Measurement residual: y = z - h(x), where h(x) is the matrix that 
    % maps the state onto the measurement.
    y = [quat(i, :) 0 0 0]' - h;
    y_output(i, :) = y';
    
    % The H matrix maps the measurement to the states.
    H=[1 0 0 0 0 0 0;
       0 1 0 0 0 0 0;
       0 0 1 0 0 0 0;
       0 0 0 1 0 0 0;
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0];
     
    % Measurement covariance update.
    S = H * P * H' + R;

    % Calculate Kalman gain.
    K = P * H' / S;

    % Corrected model prediction.
    x = x + K * y;      

    % Re-normalize Quaternion.
    qnorm = sqrt(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2);
    x(1) = x(1) / qnorm;
    x(2) = x(2) / qnorm;
    x(3) = x(3) / qnorm;
    x(4) = x(4) / qnorm;

    % Update state covariance with new knowledge.
    I = eye(length(P));
    P = (I - K * H) * P;  

    % Populate the orientation quaternion matrix.
    q(i, :) = [x(1), x(2), x(3), x(4)];
    % Populate the gyroscope bias matrix.
    wb(i, :) = [x(5), x(6), x(7)];
    
end

end
% END OF FUSIONGEKF_QUEST FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function percent = improvement_perc(nongated_rmse, gated_rmse)

% FUNCTION IMPROVEMENT_PERC computes the percentage of improvement from
% using a gating approach over a standard sensor fusion approach.
%
% - INPUT PARAMETERS:
%   |_ 'nongated_rmse': RMSE of the non gated approach.
%   |_ 'gated_rmse': RMSE of the gated approach.
%
% - OUTPUT PARAMETERS:
%   |_ 'percent': Improvement (in %).
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

x = gated_rmse * 100 / nongated_rmse;
percent = 100 - x;

end
% END OF IMPROVEMENT_PERC FUNCTION

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\