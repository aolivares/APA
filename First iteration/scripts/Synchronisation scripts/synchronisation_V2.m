% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------------ Extract Signals --------------------------------
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
% * Last modification: 25/02/2015
% -------------------------------------------------------------------------
% INFORMATION: This file contains the routine to extract, process and
% synchronise the signals from Forceplate and Gait Watch.
% Firstly, it extracts the FP signals from the ten separates forceplate 
% cycles. After, these signals are stored in a vector variable.
% After, the GaitWatch signals are processed and the most interesting of
% them are stored in a single vector (this is a modified version of the 
% 'main.m' file implemented by Alberto Olivares).
% Finally, it carries out the synchronisation between all them.
% The file is structured as follows:
% 
% * 1) Select Excel-file and obtain the GW file names.
% 
% * 2) Extract FP data from .txt files.
% 
% * 3) Load and calibrate GW data.
% 
% * 4) Synchronise signals.
% 
% * 5) Save the synchronised data in *.mat file for each person.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

% Set flags which control the visibility of the figures.
showPlotsCheck = 'yes';
showPlotsAccShank = 'yes';
showPlotsGyroShank = 'yes';
showPlotsGyroTrunk = 'no';

% Suppress warnings if no peak is detected during the calibration.
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% -------------------------------------------------------------------------
% 1) Select, read and obtain information from the excel file.
% -------------------------------------------------------------------------
 
% Select only one data files with a dialog box (only .xlsl files).
[filename_excel, filepath] = uigetfile('*.xlsx', ...
    'Select the Excel file with the data (.xlsx)', '../../data');

% Read data from *.xlsx where are stored all filenames and other 
% interesting information.
[~,file_excel] = xlsread([filepath, filename_excel]);
[rows,columns] = size(file_excel);

% The two first rows are unuseful so we obtain the data beginning from 
% second rows. We extract GW names that are in the third column, the 
% number of FP files and their names to read them afterwards.
ind = 1;

for i = 3:rows-1
        
        % Check if the cell is empty.
        if  isempty( char(strtrim(file_excel(i,4))) ) == 0  

            % We extract the GW name that identifies the patient.
            filename_GW_total (ind) = strtrim(file_excel(i,4)); 
            index_GW (ind) = i;
            ind = ind+1;
        end      
end

% Add the number of rows because it's necessary the end index to extract
% the FP signals.
index_GW = [index_GW rows+1];


% Data processing for each patient identified previously.
for j = 1:length(index_GW)-1
     
% -------------------------------------------------------------------------
% 2) Extraction FP data from  .txt files.
% -------------------------------------------------------------------------
    
% Status message.
fprintf(['Patient ', num2str(j), ' of ', num2str(length(index_GW)-1), ...
         ':', '\n\nRead data from text files...\n\n']);

% Set the firt and end index to extract the FP files names from excel.
start_data_patient = index_GW(j);
end_data_patient = index_GW(j+1)-1;

% Extract the names block.
filename_complete_FP = strtrim(file_excel(...
                               start_data_patient:end_data_patient,1));

n_files = length(filename_complete_FP);

% -------------------------------------------------------------------------
% 2.1) Read the blocks of data and compute the interesting signals.
% -------------------------------------------------------------------------

% Define the variables used to store the data.
force_sensors = cell(n_files, 1);
time_FP = cell(n_files, 1);
force_cells = cell(n_files, 1);
ML_COP = cell(n_files, 1);
AP_COP = cell(n_files, 1);
force_sum = cell(n_files, 1);
midline = cell(n_files, 1);

% Read each selected data file.
for i = 1:n_files

    % Convert the namefile in a char type to be able to read the file.
    filename_FP = char(filename_complete_FP(i));

    % Read data file.
    [header, first_line, abta, count, col_first_block] = ...
        readFPHeader(fullfile('../../data/ForcePlate/', ...
                     filename_FP));

    % Read formatted data from text file. As always "fid" is a file 
    % identifier that we obtain with "fopen".
    fid = fopen(fullfile('../../data/ForcePlate/', ...
                filename_FP));

    % Specification of the format of the data fields (double, %n) an cycles.
    % It skips the firsts lines of the data, and then reads the remaining 
    % data.
    data = textscan(fid, '%n', col_first_block * count, 'headerlines', ...
                    first_line - 1);

    % Reshape the data. This converts the data from a single column block 
    % to a five column block.
    data = reshape(data{1}, col_first_block, count)';

    % Find the starting and ending point of data. This removes all the zero
    % readings, i.e. the time instants in which the patient is not walking  
    % on the force plate.
    start_end = any(data(:, 2:5), 2);                                  
    start_end = find(diff(start_end) ~= 0)+1; 
    data_start = start_end(1);
    data_end = start_end(2);

    % Define the four pressure signals of the different platforce sensors 
    % as well as the time.
    time_FP_cycle = data(data_start:data_end, 1)';
    force_sensors_cycle = data(data_start:data_end, 2:5)';

    % Read 10 lines which contain information concerning the following data 
    % blocks.
    C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');   

    % Find out no of columns per block.
    nn = textscan(char(C2{1}(5)), '%s%n');                              
    columns = nn{2};

    % Find out no of lines per block. 
    nn = textscan(char(C2{1}(6)), '%s%n');                             
    lines = nn{2};

    % Set format to read the data, first value is a 'string' variable, the 
    % rest are 'double' variables.
    format = ['%s', repmat('%n', 1, columns)];                            

    % Read data block by block and evaluate only those which contain data
    % only to define the midline. 
    m = zeros(columns,1);           

    for b = 1:data_end

         % Read 8 lines which contain some irrelevant information.
         CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
         C3 = textscan(fid, format, lines);

         if b >= data_start
            for c = 1:columns
                % Sum up the mean of all columns.
                m(c) = m(c) + mean(C3{c + 1});                             
            end
         end
    end

    % Close data file.
    fclose(fid);

    % Definition of the midline between both feet.
    [pks, locs] = findpeaks(m, 'minpeakdistance', 10, 'minpeakheight', 100);

    if length(pks) == 2
        [~,ml] = min(m(locs(1):locs(2)));
        ml = ml + locs(1) - 1;
    else
        disp('midline could not be defined');
        beep;
    end

    % Read data one more time.
    % Define the numbers of interesting data.
    n = data_end - data_start + 1;

    % 3 values per frame(x,y,N), for 3 feet: No 1 is right foot, 2 is left 
    % foot and 3 is both feet.

    % Values: (1): COP (Center of Gravity) X (right-left).
    %         (2): COP Y (anteroposterior i.e. front-to-back).
    %         (3): Posterior margin (only Y).
    %         (4): Anterior margin (only Y is loaded).
    %         (5): Overall pressure.

    % Define the variable that contains the force data.
    force_cell = zeros(lines, columns);
    force_cell_complete_cycle = zeros(lines, columns, n);

    fid = fopen(fullfile('/APA/First iteration/data/ForcePlate/',...
                filename_FP));
    dummy_data = textscan(fid, '%n', col_first_block * count, 'headerlines',...
                 first_line - 1);
    C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');

    % Skip these data they contain no values.
    for b = 1:data_start - 1  

         % Read 8 lines which contain some irrelevant information.
         CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
         C3 = textscan(fid, format, lines);
    end

    % Define the variable that contains the ML, AP and force data.
    lateral_COP = zeros(n,3);
    antpost_COP_cycle = zeros(n,3);
    force_complete_cycle = zeros(n,3);

    for b = data_start:data_end

        % Index for data.
        c = b - data_start + 1; 

        % Read 8 lines which contain some irrelevant information.
        CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
        C3 = textscan(fid, format, lines);
        for cc = 1:columns
           force_cell(:, cc) = C3{cc + 1};
        end

        % Store the matriz with the force of each cell in the correct.
        % position.
        force_cell_complete_cycle(:,:, n) = force_cell;

        % Calculate COP in the differents cases.
        % Side 1 right, side 2 left, side 3 both.
        for side = 1:3                                                              
            switch side
                case 1
                    an = 1; 
                    en = ml;
                case 2
                    an = ml + 1;
                    en = columns;                                            
                case 3
                    an = 1; 
                    en = columns;
            end

            % Column vector. Calculate COP X to represent the markers.
            m = sum(force_cell(:, an:en));   

            % Check whether any data is >0.
            if any(m)    
                % Put the COP X in the rigth position. We multiply this  
                % result by 8.5 because each cell has 8.5 mm (as much width
                % as heigth).
                lateral_COP(c,side) = 8.5.*sum(m .* (an:en)) ./ sum(m);

                % Calculate COP Y.
                m = sum(force_cell(:, an:en), 2);                                          
                antpost_COP_cycle(c, side) = 8.5.*sum(m .* (1:lines)') ...
                                             ./sum(m);

                % Calculate rearmost foot pressure point (heel position) 
                % (y only).
                force_complete_cycle(c, side) = sum(m);
            end
        end
    end


    % Close data file.
    fclose(fid);

    % Define the 0 number like not a number.
    lateral_COP(lateral_COP == 0) = NaN;
    antpost_COP_cycle(antpost_COP_cycle == 0) = NaN;
    force_complete_cycle(force_complete_cycle == 0) = NaN;

   % Definition of the midline in mm.
    midline_mm = 8.5*ml;

   % Calculate COP X with respect to midline.
    medlateral_COP_cycle = lateral_COP - midline_mm; 


    % Store the data of each cycle.
    force_sensors{i} = force_sensors_cycle;
    time_FP{i} = time_FP_cycle./1000;
    force_cells{i} = force_cell_complete_cycle;
    ML_COP{i} = medlateral_COP_cycle';
    AP_COP{i} = antpost_COP_cycle';
    force_sum{i} = force_complete_cycle';
    midline{i} = midline_mm;

end

% Clear the unuseful variables.
clearvars -except force_sensors time_FP force_cells ML_COP ...
AP_COP force_sum midline filename_GW_total index_GW filename_FP ...
file_excel j showPlotsCheck showPlotsAccShank showPlotsGyroTrunk ...
showPlotsGyroShank


% -------------------------------------------------------------------------
% 3) Calibrate GW data.
% -------------------------------------------------------------------------
gw = gwLibrary;

fprintf('\nCalibration in progress...\n');

% Extract the GW file name for the each patient.
filename_GW = filename_GW_total(j);


% -------------------------------------------------------------------------
% 3.1) Define structure of data array.
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

size_data_struct = size(data_struct);    

% -------------------------------------------------------------------------
% 3.2) Load data file.
% -------------------------------------------------------------------------

% Load data from the hard drive.
load(fullfile('../../data/GaitWatch/',char(filename_GW)));

% Build vector containing time samples.
[f, date, start_time, end_time, file_id] = gw.getFHinfo(FileHeader);

% Build time signal. 
len_data = length(data);
time_GW = (0:len_data-1) / f;

% Reshape data (split channel 23 into 3 channels and append them to the
% data matrix).
[mag_x, mag_y, mag_z] = gw.getMagData(double(data(:,23)));

% Compute the sampling frequency of the magnetometer.
f_mag = f/(length(data)/length(mag_x));

% We now build a time signal using the newly computed frequency. This
% signal will only be used in figures (to plot magnetic field vs. time). 
time_mag=zeros(1,length(mag_x));
for i=1:length(mag_x)-1
 time_mag(i+1) = time_mag(i)+1/f_mag;
end

% -------------------------------------------------------------------------
% 3.3) Correct magnetometer signals automatically.
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

% Interpolate the magnetometer signals.
mag_x_interp = interp1(time_mag,mag_x,time_GW,'spline');
mag_y_interp = interp1(time_mag,mag_y,time_GW,'spline');
mag_z_interp = interp1(time_mag,mag_z,time_GW,'spline');

% Remove channel 23 from the data matrix and add three new channels 
% containing the interpolated magnetometer signals to the data matrix.
data = [double(data(:,1:22)) mag_x_interp' mag_y_interp' mag_z_interp'];

% -------------------------------------------------------------------------
% 3.4) Identify and extract data.
% -------------------------------------------------------------------------
% We now have to extract the data from the 'data' matrix. This is done 
% dynamically, so the code is flexible and only depends on 'data_struct'. 
% A vector is created for each one of the magnitudes which were measured.

% Define the list of segments which are calculated. 
S = cell(1,size_data_struct(1));
for i = 1:size_data_struct(1)
    S{i} = [data_struct{i,4},' ',data_struct{i,5}];
end
S = unique(S,'stable');

% Denition of all segments.
Selection = [1 2 3 4 5 6 7];

% We now extract all the data channels and dinamically generate the
% variable name according to the information in 'data_struct'.

size_data_struct = size(data_struct);

for i = 1:length(Selection)

    % We split the selection in a two-elements cell.
    position_segment = strsplit(S{Selection(i)});
    position = position_segment{1};
    segment = position_segment{2};

    for k = 1:size_data_struct(1)

        % We now want to extract any data which were measured for the
        % selected segment(s) and position(s).
        if strcmpi(data_struct{k,4},position) && ...
                strcmpi(data_struct{k,5},segment)

                % We build the name of the variable following this format:
                % 'magnitude_axis_position_segment_calibrationType'.
                var_name = strcat(data_struct{k,2},'_',data_struct{k,3},...
                         '_',data_struct{k,4},'_',data_struct{k,5},'_',...
                num2str(data_struct{k,6}));

                % We build the code that we wish to execute.
                code_line = strcat(var_name,'=data(:,%d);');

                % We evaluate the code line.
                eval(sprintf(code_line,k));
        end
    end
end

% -------------------------------------------------------------------------
% 3.5) Calibrate acceleration.
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
% 3.6) Calibrate angular rate.
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
% 3.7) Calibrate magnetic field.
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
% 3.8) Compute orientation.
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
            pitch_acc_right_shank = (gw.correct_quad_shifts(...
                                    pitch_acc_right_shank,'deg') - 90)';

            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_right_shank = gw.fusion_KF(gy, pitch_acc_right_shank,...
                                             f, var(pitch_acc_right_shank),... 
                                             var(pitch_acc_right_shank),...
                                             var(gy), alpha_KF, beta_KF,...
                                             pitch_acc_right_shank(1))';

            % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                                    threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                                                        input_signal,...
                                                        lwin_fsd,shift_fsd);

            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_right_shank = gw.fusionGKF(gy,pitch_acc_right_shank,...
                                    f,var(pitch_acc_right_shank),...
                                    var(pitch_acc_right_shank), var(gy),...
                                    alpha1,alpha2,beta1,beta2,...
                                    pitch_acc_right_shank(1), marker_fsd)';

                % Integrate angular rate (just for comparation purposes).            
                ini_pos = pitch_acc_right_shank(1);
                pitch_gyro_right_shank = gw.integRate(1/f,gy,ini_pos);
            

        case 'left shank'
            % Define the necessary signals.
            ax = a_X_left_shank_1_C;
            az = a_Z_left_shank_1_C;
            gy = g_Y_left_shank_1_C;

            % Compute pitch using acceleration.
            pitch_acc_left_shank = atan2d(az,ax);
            pitch_acc_left_shank = (gw.correct_quad_shifts(...
                                    pitch_acc_left_shank,'deg') - 90)';

            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_left_shank = gw.fusion_KF(gy, pitch_acc_left_shank,...
                                    f, var(pitch_acc_left_shank),...
                                    var(pitch_acc_left_shank), var(gy),...
                                    alpha_KF, beta_KF,...
                                    pitch_acc_left_shank(1))';

             % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+az.^2);
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,512,...
                                   threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,T_fsd,...
                                                        input_signal,...
                                                        lwin_fsd,shift_fsd);

            % Estimate pitch using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_left_shank = gw.fusionGKF(gy,pitch_acc_left_shank,...
                                    f,var(pitch_acc_left_shank),...
                                    var(pitch_acc_left_shank),...
                                    var(gy),alpha1,alpha2,beta1,beta2,...
                                    pitch_acc_left_shank(1), marker_fsd)';

                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_left_shank(1);
                pitch_gyro_left_shank = gw.integRate(1/f,gy,ini_pos);


        case 'right thigh'
            % Define the necessary signals.
            ax = a_X_right_thigh_1_C;
            az = a_Z_right_thigh_1_C;
            gy = g_Y_right_thigh_1_C;

            % Compute pitch using acceleration.
            pitch_acc_right_thigh = atan2d(az,ax);
            pitch_acc_right_thigh = (gw.correct_quad_shifts(...
                                     pitch_acc_right_thigh,'deg') - 90)';

            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_right_thigh = gw.fusion_KF(gy, pitch_acc_right_thigh,...
                                     f, var(pitch_acc_right_thigh), ...
                                     var(pitch_acc_right_thigh),var(gy), ...
                                     alpha_KF, beta_KF,...
                                     pitch_acc_right_thigh(1))';

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
                                    f,var(pitch_acc_right_thigh),...
                                    var(pitch_acc_right_thigh),var(gy),...
                                    alpha1,alpha2,beta1,beta2,...
                                    pitch_acc_right_thigh(1), marker_fsd)';

                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_right_thigh(1);
                pitch_gyro_right_thigh = gw.integRate(1/f,gy,ini_pos);

        case 'left thigh'
            % Define the necessary signals.
            ax = a_X_left_thigh_1_C;
            az = a_Z_left_thigh_1_C;
            gy = g_Y_left_thigh_1_C;

            % Compute pitch using acceleration.
            pitch_acc_left_thigh = atan2d(az,ax);
            pitch_acc_left_thigh = gw.correct_quad_shifts(...
                                    pitch_acc_left_thigh,'deg');
            pitch_acc_left_thigh = (pitch_acc_left_thigh - 90)';

            % Fuse acceleration-based pitch with angular rate using 
            % Kalman filter.
            alpha_KF = 1000;
            beta_KF = 0.001;
            pitch_KF_left_thigh = gw.fusion_KF(gy, pitch_acc_left_thigh,...
                                  f, var(pitch_acc_left_thigh) ,...
                                  var(pitch_acc_left_thigh), var(gy), ...
                                  alpha_KF, beta_KF,...
                                  pitch_acc_left_thigh(1))';

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
                                   f,var(pitch_acc_left_thigh),...
                                   var(pitch_acc_left_thigh),...
                                   var(gy),alpha1,alpha2,beta1,beta2,...
                                   pitch_acc_left_thigh(1), marker_fsd)';

                % Integrate angular rate (just for comparation purposes).
                ini_pos = pitch_acc_left_thigh(1);
                pitch_gyro_left_thigh = gw.integRate(1/f,gy,ini_pos);

        case 'left arm'
            % Define the necessary signals.
            gx = g_X_left_arm_1_C;
            gy = g_Y_left_arm_1_C;

            % Integrate angular rate to compute pitch and roll.
            pitch_gyro_left_arm = gw.integRate(1/f,gy,0);
            roll_gyro_left_arm = gw.integRate(1/f,gx,0);

            % Apply high pass filter to partially remove bias.
            lower_freq_limit=0.2; 
            [b,a]=butter(3,lower_freq_limit*2/(f),'high'); 
            pitch_gyro_left_arm_f = (filtfilt(b,a,pitch_gyro_left_arm) + ...
                                     pitch_gyro_left_arm(1));
            roll_gyro_left_arm_f  = (filtfilt(b,a,roll_gyro_left_arm) + ...
                                    roll_gyro_left_arm(1));



        case 'right arm'
            % Define the necessary signals.
            gx = g_X_right_arm_1_C;
            gy = g_Y_right_arm_1_C;

            % Integrate angular rate to compute pitch and roll.
            pitch_gyro_right_arm = gw.integRate(1/f,gy,0);
            roll_gyro_right_arm = gw.integRate(1/f,gx,0); 

            % Apply high pass filter to partially remove bias.
            lower_freq_limit=0.2; 
            [b,a] = butter(3,lower_freq_limit*2/(f),'high'); 
            pitch_gyro_right_arm_f = filtfilt(b,a,pitch_gyro_right_arm)+...
                                     pitch_gyro_right_arm(1);
            roll_gyro_right_arm_f  = filtfilt(b,a,roll_gyro_right_arm)+...
                                     roll_gyro_right_arm(1);


        case 'center trunk'
            % Define the necessary signals.
            delay = 1;
            ax = a_X_center_trunk_3_C(delay:end);
            ay = a_Y_center_trunk_3_C(delay:end);
            az = a_Z_center_trunk_3_C(delay:end);
            gx = g_X_center_trunk_1_C(delay:end);
            gy = g_Y_center_trunk_1_C(delay:end);
            gz = g_Z_center_trunk_1_C(delay:end);
            hx = h_X_center_trunk_3_C(delay:end);
            hy = h_Y_center_trunk_3_C(delay:end);
            hz = h_Z_center_trunk_3_C(delay:end);  

            % Compute roll, pitch using accelerometer.
            pitch_acc = atan2d(sqrt(ay.^2+az.^2),ax) - 90;
            roll_acc = atan2d(ay,az);

            % Find the projections of hx and hy in the XY plane.
            hx_n = -hx.*cosd(pitch_acc) + hy.*sind(pitch_acc).*...
                    sind(roll_acc) - hz.*sind(pitch_acc).*cosd(roll_acc);
            hy_n = hy.*cosd(roll_acc) + hz.*sind(roll_acc);

            % Compute yaw.
            yaw_mag = atan2d(hx_n,hy_n);

            % Compensate quadrante discontinuities.
            yaw_mag = gw.correct_yaw_quad_shifts(yaw_mag,'deg');

            % Compute roll and pitch using non-quaternion Kalman Filter. 
            alpha_KF = 100;
            beta_KF = 0.01;
            pitch_KF_center_trunk = gw.fusion_KF(gy,pitch_acc,f,...
                                    var(pitch_acc),var(pitch_acc),var(gy),...
                                    alpha_KF,beta_KF,pitch_acc(1))';
            roll_KF_center_trunk = gw.fusion_KF(gx,roll_acc,f,...
                                    var(roll_acc),var(roll_acc),var(gx),...
                                    alpha_KF,beta_KF, roll_acc(1))';
            alpha_KF = 1000;
            beta_KF = 0.001;
            yaw_KF_center_trunk = gw.fusion_KF(gz,yaw_mag,f,...
                                     var(yaw_mag),var(yaw_mag),var(gz),...
                                     alpha_KF,beta_KF, yaw_mag(1))';


            % Find the projections of hx and hy in the XY plane using 
            % fused pitch and roll.
            hx_n_KF = -hx.*cosd(pitch_KF_center_trunk) + ...
                      hy.*sind(pitch_KF_center_trunk).*...
                      sind(roll_KF_center_trunk) - hz.*...
                      sind(pitch_KF_center_trunk).*cosd(roll_KF_center_trunk);
            hy_n_KF = hy.*cosd(roll_KF_center_trunk) +...
                      hz.*sind(roll_KF_center_trunk);

            % Compute yaw.
            yaw_mag_kf = atan2d(hx_n_KF,hy_n_KF);

            % Compensate quadrante discontinuities.
            yaw_mag_kf = gw.correct_yaw_quad_shifts(yaw_mag_kf,'deg');
            correct_again = 1;
            limit = 250;
            while correct_again == 1
                correct_again = 0;
                for l = 2: length(yaw_mag_kf)
                    if abs(yaw_mag_kf(l)-yaw_mag_kf(l-1)) > limit
                        correct_again = 1;
                    end
                end
                if correct_again == 1
                    yaw_mag_kf = gw.correct_yaw_quad_shifts(...
                                 yaw_mag_kf,'deg')';
                end
            end          
            yaw_2KF_center_trunk = gw.fusion_KF(gz,yaw_mag_kf,f,...
                                   var(yaw_mag_kf),var(yaw_mag_kf),var(gz),...
                                   alpha_KF, beta_KF,yaw_mag_kf(1))';

            % Compute intensity level.         
            lwin_fsd = 20;    
            threshold_fsd = 3;    
            shift_fsd = 19;    
            lambda = 30;
            input_signal = sqrt(ax.^2+ay.^2+az.^2)';
            [V_fsd,T_fsd] = gw.fsd(input_signal,lwin_fsd,shift_fsd,...
                            512,threshold_fsd);
            [marker_fsd,T_fsd_expanded] = gw.compEstMark(V_fsd,...
                                           T_fsd,input_signal,lwin_fsd,...
                                           shift_fsd);

            % Estimate Euler angles using Gated Kalman filter.
            alpha1 = 100;
            alpha2 = 10000;
            beta1 = 0.001;
            beta2 = 0.00001;
            pitch_GKF_center_trunk = gw.fusionGKF(gy,pitch_acc,f,...
                                     var(pitch_acc),var(pitch_acc),...
                                     var(gy),alpha1,alpha2, beta1,beta2,...
                                     pitch_acc(1),marker_fsd)';
            roll_GKF_center_trunk = gw.fusionGKF(gx,roll_acc,f,...
                                    var(roll_acc),var(roll_acc),var(gx),...
                                    alpha1,alpha2, beta1,beta2,roll_acc(1),...
                                    marker_fsd)';
            yaw_GKF_center_trunk = gw.fusionGKF(gz,yaw_mag,f,...
                                   var(yaw_mag),var(yaw_mag),var(gz),...
                                   alpha1,alpha2, beta1,beta2,yaw_mag(1),...
                                   marker_fsd)';
            yaw_2GKF_center_trunk = gw.fusionGKF(gz,yaw_mag_kf,f,...
                                    var(yaw_mag_kf), var(yaw_mag_kf),...
                                    var(gz),alpha1,alpha2,beta1,beta2,...
                                    yaw_mag_kf(1),marker_fsd)';

            % Compute roll, pitch and yaw using quaternion Extended 
            % Kalman Filter.
            gyroVarX = 0.1;
            gyroVarY = 0.1;
            gyroVarZ = 0.1;
            mu_gain = 5;
            alpha = 5;
            q_ini = gw.eulerToQuat(roll_acc(1)/180*pi,pitch_acc(1)/...
                                    180*pi, yaw_mag(1)/180*pi);

            state_ini = [q_ini(1), q_ini(2), q_ini(3), 0.5]';

            [roll_EKF_center_trunk, pitch_EKF_center_trunk, ...
             yaw_EKF_center_trunk] = gw.quat9dofEKF(ax, ay, az,...
                                        gx, gy, gz, hx, hy, -hz, gyroVarX,...
                                        gyroVarY, gyroVarZ, alpha, mu_gain,...
                                        f, state_ini);           

                pitch_gyro = gw.integRate(1/f,gy,pitch_acc(1));
                roll_gyro = gw.integRate(1/f,gx,roll_acc(1)); 
                yaw_gyro = gw.integRate(1/f,gz,yaw_mag(1)); 
      end
end

% Turn some signals into row vectors to use them easier afterwards.
a_X_left_shank_1_C =a_X_left_shank_1_C';
a_X_left_thigh_1_C =a_X_left_thigh_1_C';
a_X_right_shank_1_C =a_X_right_shank_1_C';
a_X_right_thigh_1_C =a_X_right_thigh_1_C';

a_Z_left_shank_1_C =a_Z_left_shank_1_C';
a_Z_left_thigh_1_C =a_Z_left_thigh_1_C';
a_Z_right_shank_1_C =a_Z_right_shank_1_C';
a_Z_right_thigh_1_C=a_Z_right_thigh_1_C';

% Subtract the mode of the gyroscope's signals to centrate these signals in
% the 0 value.
g_X_center_trunk_1_C = g_X_center_trunk_1_C' - mode(g_X_center_trunk_1_C);
g_X_left_arm_1_C = g_X_left_arm_1_C' - mode(g_X_left_arm_1_C);
g_X_right_arm_1_C = g_X_right_arm_1_C' - mode(g_X_right_arm_1_C);
g_Y_center_trunk_1_C = g_Y_center_trunk_1_C' - mode(g_Y_center_trunk_1_C);
g_Y_left_shank_1_C= g_Y_left_shank_1_C' - mode(g_Y_left_shank_1_C);
g_Y_left_arm_1_C = g_Y_left_arm_1_C' - mode(g_Y_left_arm_1_C);
g_Y_left_thigh_1_C = g_Y_left_thigh_1_C' - mode(g_Y_left_thigh_1_C);
g_Y_right_shank_1_C = g_Y_right_shank_1_C' - mode (g_Y_right_shank_1_C);
g_Y_right_arm_1_C = g_Y_right_arm_1_C' - mode(g_Y_right_arm_1_C);
g_Y_right_thigh_1_C = g_Y_right_thigh_1_C' - mode(g_Y_right_thigh_1_C);
g_Z_center_trunk_1_C = g_Z_center_trunk_1_C' - mode(g_Z_center_trunk_1_C);  

% ---------------------------------------------------------------------
% 3.9) Clear variables
% ---------------------------------------------------------------------

clearvars -except force_sensors time_FP force_cells ML_COP ...
AP_COP force_sum midline filename_GW_total filename_FP index_GW ...
file_excel j showPlotsCheck showPlotsAccShank showPlotsGyroTrunk ...
showPlotsGyroShank ...
a_X_center_trunk_3_C a_X_left_shank_1_C a_X_left_thigh_1_C...
a_X_right_shank_1_C a_X_right_thigh_1_C a_Y_center_trunk_3_C ...
a_Z_center_trunk_3_C a_Z_left_shank_1_C a_Z_left_thigh_1_C...
a_Z_right_shank_1_C a_Z_right_thigh_1_C...
g_X_center_trunk_1_C g_X_left_arm_1_C g_X_right_arm_1_C ...
g_Y_center_trunk_1_C g_Y_left_shank_1_C g_Y_left_arm_1_C ...
g_Y_left_thigh_1_C g_Y_right_shank_1_C g_Y_right_arm_1_C ...
g_Y_right_thigh_1_C g_Z_center_trunk_1_C...
h_X_center_trunk_3_C h_Y_center_trunk_3_C h_Y_center_trunk_3_C...
pitch_acc_left_shank pitch_gyro_left_shank pitch_KF_left_shank...
pitch_GKF_left_shank ...
pitch_acc_right_shank pitch_gyro_right_shank ...
pitch_KF_right_shank pitch_GKF_right_shank ...
pitch_acc_left_thigh pitch_gyro_left_thigh pitch_KF_left_thigh ...
pitch_GKF_left_thigh pitch_GKF_right_thigh...
pitch_acc_right_thigh pitch_gyro_right_thigh ...
pitch_KF_right_thigh pitch_GKF_left_thigh ...
pitch_gyro_left_arm pitch_gyro_left_arm_f pitch_gyro_right_arm ...
pitch_gyro_right_arm_f roll_gyro_left_arm roll_gyro_left_arm_f ...
roll_gyro_right_arm roll_gyro_right_arm_f ...
pitch_acc pitch_gyro pitch_KF_center_trunk ...
pitch_GKF_center_trunk ...
roll_acc roll_gyro roll_KF_center_trunk ...
roll_GKF_center_trunk ...
yaw_mag_k yaw_gyro yaw_KF_center_trunk yaw_GKF_center_trunk ...
yaw_2KF_center_trunk yaw_2GKF_center_trunk time_GW
      
 
% -------------------------------------------------------------------------
% 4) Synchronisation.
% -------------------------------------------------------------------------

% Status message.
fprintf('\nSynchronisation in progress...\n');

% Set sampling frequencies of the Force Plate and the GaitWatch.
fs_FP = 120;
fs_GW = 200;

% -------------------------------------------------------------------------
% 4.1) Application of Activity detectors
% -------------------------------------------------------------------------
wag = wagLibrary;    

% Define input signal.
axC = a_X_right_shank_1_C;
azC = a_Z_right_shank_1_C;

input_signal = sqrt(axC .^ 2 + azC .^ 2)';

% Computation of intensity markers.

% FSD (window size, decision threshold, overlapping and normalization
% factor). 
lwin_fsd = 100;  threshold_fsd = 3;  shift_fsd = 100; lambda = 50;

% LTSD (window size, decision threshold and overlapping).
lwin_ltsd = 100;       threshold_ltsd = 4.3;   shift_ltsd = 10;

% 3) Get the decision signal of the FSD algorithm and the marker.
[V_fsd, T_fsd] = wag.fsd(input_signal, lwin_fsd, shift_fsd, 512, ...
    threshold_fsd);
[marker_fsd, T_fsd_expanded] = wag.compEstMark(V_fsd, T_fsd, input_signal, ...
    lwin_fsd, shift_fsd);

% 4) Get the decision signal of the LTSD algorithm and the marker.
[V_ltsd, T_ltsd] = wag.ltsd(input_signal, lwin_ltsd, shift_ltsd, 512, ...
    threshold_ltsd);
[marker_ltsd, T_ltsd_expanded] = wag.compEstMark(V_ltsd, T_ltsd, ...
    input_signal, lwin_ltsd, shift_ltsd);

if strcmpi(showPlotsCheck,'yes')    
figure
subplot(2, 1, 1)
plot(T_fsd_expanded)
hold on
plot(threshold_fsd * ones(1, length(T_fsd_expanded)), 'r')
legend('Detector output (FSD)', 'Detection threshold')
subplot(2, 1, 2)
plot(T_ltsd_expanded)
hold on
plot(threshold_ltsd * ones(1, length(T_ltsd_expanded)), 'r')
legend('Detector output (LTSD)', 'Detection threshold')

figure
subplot(2, 1, 1)
plot(input_signal)
hold on
plot(marker_fsd + 1, 'r')
legend('Input signal','FSD decision')
subplot(2, 1, 2)
plot(input_signal)
hold on
plot(marker_ltsd + 1, 'r')
legend('Input signal','LTSD decision')
end
% -------------------------------------------------------------------------
% 4.1) Find the second positive peak of each cycle in the GaitWatch signal 
%      of the z-axis of the acceleration of the shank, that is, the point 
%      in time when the patient walks onto the force plate. 
% -------------------------------------------------------------------------

% Set tuning parameter for peak detection in acceleration signal.
threshold_pos = 1.05;    
threshold_neg = -0.9;

% Determinate the initial and end point of each interval where we need to 
% find the peaks, i.e, the first activity period of each cycle. 
edges = find(diff(marker_ltsd)~=0);
initcross = edges(1:4:length(edges));
finalcross = edges(2:4:length(edges));

% Calculate the peaks in each interval.
for k = 1:length(initcross)

    % Find all peaks smaller than threshold in each interval in the left 
    % shank signal.
    [neg_peak_values_l, neg_peak_locations_l] = findpeaks(...
                            -a_Z_left_shank_1_C(initcross(k):finalcross(k)),...
                            'minpeakheight', threshold_neg);

    % Store the index of the highest negative peak.                                      
    neg_peaks_l(k) = find(a_Z_left_shank_1_C(initcross(k):finalcross(k))== ...
                    -max(neg_peak_values_l), 1, 'last') + initcross(k) - 1;
                                     
 
    % Find all peaks greater than threshold in the interval specified by the 
    % minimum peak calculated above and the end of the interval.
    [peak_values_l, peak_locations_l] = findpeaks(a_Z_left_shank_1_C(...
                                        neg_peaks_l(k):finalcross(k)), ...
                                        'minpeakheight', threshold_pos);

    % Store the index of the highest positive peak .                                      
    sync_peaks_l(k) = find(...
                     a_Z_left_shank_1_C(neg_peaks_l(k):finalcross(k)) ...
                     == max(peak_values_l), 1, 'first') + neg_peaks_l(k) - 1;
                             
    % Find all peaks smaller than threshold in each interval in the right 
    % shank signal.
    [neg_peak_values_r, neg_peak_locations_r] = findpeaks(...
                            -a_Z_right_shank_1_C(initcross(k):finalcross(k)),...
                            'minpeakheight', threshold_neg);

    % Store the index of the highest negative peak.                                      
    neg_peaks_r(k) = find(a_Z_right_shank_1_C(initcross(k):finalcross(k))== ...
                    -max(neg_peak_values_r), 1, 'last') + initcross(k) - 1;
                                     
 
    % Find all peaks greater than threshold in the interval specified by the 
    % minimum peak calculated above and the end of the interval.
    [peak_values_r, peak_locations_r] = findpeaks(a_Z_right_shank_1_C(...
                                        neg_peaks_r(k):finalcross(k)), ...
                                        'minpeakheight', threshold_pos);

    % Store the index of the highest positive peak .                                      
    sync_peaks_r(k) = find(...
                     a_Z_right_shank_1_C(neg_peaks_r(k):finalcross(k)) ...
                     == max(peak_values_r), 1, 'first') + neg_peaks_r(k) - 1;

end

% Evaluate if the patient steps with the left or right limb first for each
% cycle and store the first peak in sync_peaks, respectively, then sort it.
% In this case, we use the acceleration data.
sync_peaks_acc = [sync_peaks_r(sync_peaks_l > sync_peaks_r), ...
              sync_peaks_l(sync_peaks_r > sync_peaks_l)];
sync_peaks_acc = sort(sync_peaks_acc);


% Store last peak of each cycle.
 last_peaks = edges(4:4:length(edges));

% -------------------------------------------------------------------------
% 4.2 Store seperate GaitWatch cycles in a cell array of time series
%     objects and add the point in time when the patient walks on the
%     force plate as an event.
% -------------------------------------------------------------------------

% Set the additional time in seconds before (sync_peaks) and after
% (last_peaks) each cycle that will be stored in the time series.
add_time = 2;
add_samples = floor(add_time * fs_GW);

% Create cell array containing the time series of the seperate cycles of 
% the acceleration of the trunk.
a_trunk = createTimeseriesGW([a_X_center_trunk_3_C; a_Y_center_trunk_3_C; ...
                              a_Z_center_trunk_3_C], time_GW, sync_peaks_acc, ...
                              last_peaks, add_samples, ...
                              'Acceleration trunk', 'seconds', 'g');
                        
% Create cell array containing the time series of the seperate cycles of 
% the acceleration of the thighs.
a_thighs = createTimeseriesGW([a_X_left_thigh_1_C;  a_Z_left_thigh_1_C; ...
                               a_X_right_thigh_1_C; a_Z_right_thigh_1_C], ...
                               time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                               'Acceleration thighs', 'seconds', 'g');
                        
% Create cell array containing the time series of the seperate cycles of 
% the acceleration of the shanks.
a_shanks = createTimeseriesGW([a_X_left_shank_1_C;  a_Z_left_shank_1_C; ...
                               a_X_right_shank_1_C; a_Z_right_shank_1_C], ...
                               time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                               'Acceleration shanks', 'seconds', 'g');
                          
                          
% Create cell array containing the time series of the seperate cycles of 
% the angular rate of the trunk.
g_trunk = createTimeseriesGW([g_X_center_trunk_1_C; g_Y_center_trunk_1_C; ...
                              g_Z_center_trunk_1_C], time_GW, sync_peaks_acc, ...
                              last_peaks, add_samples, ...
                              'Angular rate trunks', 'seconds', 'deg/s');
                             
% Create cell array containing the time series of the seperate cycles of 
% the angular rate of the thighs.
g_thighs = createTimeseriesGW([g_Y_left_thigh_1_C; g_Y_right_thigh_1_C], ...
                               time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                               'Angular rate thighs', 'seconds', 'deg/s');
                             
% Create cell array containing the time series of the seperate cycles of 
% the angular rate of the shanks.
g_shanks = createTimeseriesGW([g_Y_left_shank_1_C; g_Y_right_shank_1_C], ...
                               time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                               'Angular rate shanks', 'seconds', 'deg/s');
                             
% Create cell array containing the time series of the seperate cycles of 
% the angular rate of the arms.
g_arms = createTimeseriesGW([g_X_left_arm_1_C; g_Y_left_arm_1_C; ...
                             g_X_right_arm_1_C; g_Y_right_arm_1_C], ...
                             time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                             'Angular rate arms', 'seconds', 'deg/s');
                             
                             
% Create cell array containing the time series of the seperate cycles of 
% the magnetic field at the trunk.
h_trunk = createTimeseriesGW([h_X_center_trunk_3_C; h_Y_center_trunk_3_C], ...
                              time_GW, sync_peaks_acc, last_peaks, add_samples, ...
                              'Magnetic field trunk', 'seconds', 'Gauss');
                          

% Create cell array containing the time series of the seperate cycles of 
% the pitch of the thighs (computed from acceleration).
pitch_thighs_acc = createTimeseriesGW([pitch_acc_left_thigh; ...
                                       pitch_acc_right_thigh], ...
                                       time_GW, sync_peaks_acc, last_peaks, ...
                                       add_samples, ['Pitch thighs ', ...
                                       '(computed from acceleration)'], ...
                                       'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the shanks (computed from acceleration).
pitch_shanks_acc = createTimeseriesGW([pitch_acc_left_shank; ...
                                       pitch_acc_right_shank], ...
                                       time_GW, sync_peaks_acc, last_peaks, ...
                                       add_samples, ['Pitch shanks ', ...
                                       '(computed from acceleration)'], ...
                                       'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the trunk (GKF).
pitch_trunk_GKF = createTimeseriesGW(pitch_GKF_center_trunk, ...
                                     time_GW, sync_peaks_acc, last_peaks, ...
                                     add_samples, 'Pitch trunk (GKF)', ...
                                     'seconds', 'deg');

% Create cell array containing the time series of the seperate cycles of 
% the pitch of the thighs.
pitch_thighs_GKF = createTimeseriesGW([pitch_GKF_left_thigh; ...
                                       pitch_GKF_right_thigh], ...
                                      time_GW, sync_peaks_acc, last_peaks, ...
                                      add_samples, 'Pitch thighs (GKF)', ...
                                      'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the shanks (GKF).
pitch_shanks_GKF = createTimeseriesGW([pitch_GKF_left_shank; ...
                                       pitch_GKF_right_shank], ...
                                       time_GW, sync_peaks_acc, last_peaks, ...
                                       add_samples, 'Pitch thighs (GKF)', ...
                                       'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the arms (computed from gyroscope).
pitch_arms_gyro = createTimeseriesGW([pitch_gyro_left_arm; ...
                                      pitch_gyro_right_arm], ...
                                      time_GW, sync_peaks_acc, last_peaks, ...
                                      add_samples, ['Pitch arms ', ...
                                      '(computed from gyroscope)'], ...
                                      'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the thighs (computed from gyroscope).
pitch_thighs_gyro = createTimeseriesGW([pitch_gyro_left_thigh; ...
                                        pitch_gyro_right_thigh], ...
                                        time_GW, sync_peaks_acc, last_peaks, ...
                                        add_samples, ['Pitch thighs ', ...
                                        '(computed from gyroscope)'], ...
                                        'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the shanks (computed from gyroscope).
pitch_shanks_gyro = createTimeseriesGW([pitch_gyro_left_shank; ...
                                        pitch_gyro_right_shank], ...
                                        time_GW, sync_peaks_acc, last_peaks, ...
                                        add_samples, ['Pitch shanks ', ...
                                        '(computed from gyroscope)'], ...
                                        'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the trunk (KF).
pitch_trunk_KF = createTimeseriesGW(pitch_KF_center_trunk, time_GW, ...
                                    sync_peaks_acc, last_peaks, add_samples, ...
                                    'Pitch trunk (KF)', 'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the thighs (KF).
pitch_thighs_KF = createTimeseriesGW([pitch_KF_left_thigh; ...
                                      pitch_KF_right_thigh], ...
                                      time_GW, sync_peaks_acc, last_peaks, ...
                                      add_samples, 'Pitch thighs (KF)', ...
                                      'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the pitch of the shanks (KF).
pitch_shanks_KF = createTimeseriesGW([pitch_KF_left_shank; ...
                                      pitch_KF_right_shank], ...
                                      time_GW, sync_peaks_acc, last_peaks, ...
                                      add_samples, 'Pitch shanks (KF)', ...
                                      'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the roll of the trunk (GKF).
roll_trunk_GKF = createTimeseriesGW(roll_GKF_center_trunk, time_GW, ...
                                    sync_peaks_acc, last_peaks, add_samples, ...
                                   'Roll trunk (GKF)', 'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the roll of the trunk (KF).
roll_trunk_KF = createTimeseriesGW(roll_KF_center_trunk, time_GW, ...
                                   sync_peaks_acc, last_peaks, add_samples, ...
                                   'Roll trunk (KF)', 'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the roll of the arms (GKF).
roll_arms_GKF = createTimeseriesGW([roll_gyro_left_arm; ...
                                    roll_gyro_right_arm], time_GW, ...
                                    sync_peaks_acc, last_peaks, add_samples, ...
                                    'Roll arms (gyroscope)', 'seconds', 'deg');
                         
% Create cell array containing the time series of the seperate cycles of 
% the yaw of the trunk (GKF).
yaw_trunk_GKF = createTimeseriesGW(yaw_GKF_center_trunk, time_GW, ...
                                   sync_peaks_acc, last_peaks, add_samples, ...
                                   'Yaw trunk (GKF)', 'seconds', 'deg');

% Create cell array containing the time series of the seperate cycles of 
% the yaw of the trunk (2GKF).
yaw_trunk_2GKF = createTimeseriesGW(yaw_2GKF_center_trunk, time_GW, ...
                                    sync_peaks_acc, last_peaks, add_samples, ...
                                    'Yaw trunk (2GKF)', 'seconds', 'deg');

% Create cell array containing the time series of the seperate cycles of 
% the yaw of the trunk (KF).
yaw_trunk_KF = createTimeseriesGW(yaw_KF_center_trunk, time_GW, ...
                                  sync_peaks_acc, last_peaks, add_samples, ...
                                  'Yaw trunk (KF)', 'seconds', 'deg');

% Create cell array containing the time series of the seperate cycles of 
% the yaw of the trunk (2KF).
yaw_trunk_2KF = createTimeseriesGW(yaw_2KF_center_trunk, time_GW, ...
                                   sync_peaks_acc, last_peaks, add_samples, ...
                                  'Yaw trunk (2KF)', 'seconds', 'deg');


% -------------------------------------------------------------------------
% 4.3) Synchronise and resample the seperate force plate cycles, then store 
%      them in a cell array of time series objects and add the point in
%      time when the patient walks on the force plate as an event.
% -------------------------------------------------------------------------

% Compute the points in time where the synchronisation peaks appear.
sync_peak_times = time_GW(sync_peaks_acc);

% Create cell array containing the time series of the seperate cycles of
% the four force sensors.
force_sensors_ts = createTimeseriesFP(force_sensors, time_FP, sync_peak_times, ...
                                      fs_GW, 'Force Sensors', 'seconds', 'N');
                               
% Create cell array containing the time series of the seperate cycles of
% the sum of forces.                               
force_sum_ts = createTimeseriesFP(force_sum, time_FP, sync_peak_times, ...
                                  fs_GW, 'Force sum', 'seconds', 'N');
                              
% Create cell array containing the time series of the seperate cycles of
% the force cells.                               
 force_cells_ts = createTimeseriesFP(force_cells, time_FP, sync_peak_times, ...
                                     fs_GW, 'Force cells', 'seconds', 'N');                              
                               
% Create cell array containing the time series of the seperate cycles of
% the anterior-posterior center of pressure.
AP_COP_ts = createTimeseriesFP(AP_COP, time_FP, sync_peak_times, ...
                               fs_GW, 'AP-COP', 'seconds', 'mm');
                               
% Create cell array containing the time series of the seperate cycles of
% the medio-lateral center of pressure.
ML_COP_ts = createTimeseriesFP(ML_COP, time_FP, sync_peak_times, ...
                               fs_GW, 'ML-COP', 'seconds', 'mm');


% -------------------------------------------------------------------------
% 5) Store cell arrays of time series objects as .mat file
% -------------------------------------------------------------------------

% Generate the output file name.
 name_file = textscan(filename_FP,'%s','Delimiter',',');
 name_file = char(name_file{1}{1});
 name_file = strcat( name_file,'_synchronised.mat');
%  
%  save(['../../data/Synchronised/' ...
%        name_file], 'a_shanks', 'a_thighs', 'a_trunk', 'AP_COP_ts', ...
%        'force_cells_ts', 'force_sensors_ts', 'force_sum_ts', 'g_arms', ...
%        'g_shanks', 'g_thighs', 'g_trunk', 'h_trunk', 'ML_COP_ts', ...
%        'pitch_arms_gyro', 'pitch_shanks_acc', 'pitch_shanks_GKF', ...
%        'pitch_shanks_gyro', 'pitch_shanks_KF', 'pitch_thighs_acc', ...
%        'pitch_thighs_GKF', 'pitch_thighs_gyro', 'pitch_thighs_KF', ...
%        'pitch_trunk_GKF', 'pitch_trunk_KF', 'roll_arms_GKF', ...
%        'roll_trunk_GKF', 'roll_trunk_KF', 'yaw_trunk_2GKF', ...
%        'yaw_trunk_2KF', 'yaw_trunk_GKF', 'yaw_trunk_KF');
  
  fprintf('\nSaved synchronised signals!\n\n\n');
  
% Plots.

% Detected sync-peaks.
if strcmpi(showPlotsAccShank,'yes')
figure()
subplot(2, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l), a_Z_left_shank_1_C(sync_peaks_l), 'r.');

title('Acceleration left shank with detected sync-peaks left');
xlabel('Time in s');
ylabel('Acceleration in g');

subplot(2, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r), a_Z_right_shank_1_C(sync_peaks_r), 'm.');

title('Acceleration right shank with detected sync-peaks right');   
xlabel('Time in s');
ylabel('Acceleration in g');


% Left-right detection.
figure();
subplot(2, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)), 'r.', 'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_acc), 'LineStyle',':', 'LineWidth', 2 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

title(['Acceleration of z-axis of left shank with synchronisation lines and red marker when ' ...
       'patient steps with the left foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([125, 165, -0.1, 2.1]);

subplot(2, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)), 'm.', 'markersize', 20);
 
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_acc), 'LineStyle',':', 'LineWidth', 2, 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

title(['Acceleration of z-axis of right shank with synchronisation lines and magenta marker ' ...
       'when patient steps with the right foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([125, 165, -0.1, 2.3]);



% Append all separate time series of the four force sensor signals and 
% extract data from timeseries for plot.
force_sensors_complete_ts = append(force_sensors_ts{1, :});
fs_data = force_sensors_complete_ts.data;


figure();
subplot(3, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)),...
     'r.', 'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the left shank with synchronisation lines and red marker when the ' ...
       'patient steps with the left foot first.']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([0, 300, -0.1, 2.1]);

subplot(3, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     'm.', 'markersize', 20);
 
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the right shank with synchronisation lines and magenta marker ' ...
       'when the patient steps with the right foot first.']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([0, 300, -0.1, 2.3]);

subplot(3, 1, 3)
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('Exemplary the course of the synchronised signals of the four force sensors.');
xlabel('Time in s');
ylabel('Force in N');
legend(['Left front foot', 'Left back foot', 'Right front foot', ...
    'Right back foot']);

axis([0, 300, 0, 800]);
end
if strcmpi(showPlotsCheck,'yes')
figure();
subplot(3, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)),...
     'r.', 'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', 'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the left shank with synchronisation lines and red marker when the ' ...
       'patient steps with the left foot first.']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([125, 165, -0.1, 2.1]);

subplot(3, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)),...
     'm.', 'markersize', 20);
 
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the right shank with synchronisation lines and red marker when the ' ...
       'patient steps with the left foot first.']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([125, 165, -0.1, 2.3]);

subplot(3, 1, 3)
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('Exemplary the course of the synchronised signals of the four force sensors.');
xlabel('Time in s');
ylabel('Force in N');
legend('Left fore foot', 'Left hind foot', 'Right fore foot', 'Right hind foot');

axis([125, 165, 0, 1000]);

figure();
subplot(3, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)), 'r.',...
     'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the left shank with synchronisation lines and red marker when the' ...
       'patient steps with the left foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([127, 129, -0.1, 2.1]);

subplot(3, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     'm.', 'markersize', 20);
 
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the right shank with synchronisation lines and magenta marker ' ...
       'when the patient steps with the right foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([127, 129, -0.1, 2.3]);

subplot(3, 1, 3)
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('Exemplary the synchronised trace of the force sensor of the left front foot');
xlabel('Time in s');
ylabel('Force in N');

axis([127, 129, 0, 800]);


figure();
subplot(3, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l(sync_peaks_r > sync_peaks_l)), ...
     a_Z_left_shank_1_C(sync_peaks_l(sync_peaks_r > sync_peaks_l)), 'r.', ...
     'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the left shank with synchronisation lines and red marker when the' ...
       'patient steps with the left foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([150, 152, -0.1, 2.1]);

subplot(3, 1, 2)
plot(time_GW, a_Z_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r(sync_peaks_l > sync_peaks_r)), ...
     a_Z_right_shank_1_C(sync_peaks_r(sync_peaks_l > sync_peaks_r)),...
     'm.', 'markersize', 20);
 
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the z-axis of the right shank with synchronisation lines and magenta marker ' ...
       'when the patient steps with the right foot first']);
xlabel('Time in s');
ylabel('Acceleration in g');
axis([150, 152, -0.1, 2.3]);

subplot(3, 1, 3)
plot(force_sensors_complete_ts.time, reshape(fs_data(:, 1, :), ...
     [4, max(size(fs_data))]));


% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_l), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(sync_peaks_r), 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title('Exemplary the synchronised trace of the force sensor of the left front foot');
xlabel('Time in s');
ylabel('Force in N');

axis([150, 152, 0, 800]);
end

% Synchronisation with the gyroscope signal.
input_signal = sqrt(g_Y_right_shank_1_C .^ 2 + g_Y_right_shank_1_C .^ 2)';
%input_signal= g_Y_right_shank_1_C';

% FSD (window size, decision threshold, overlapping and normalization
% factor). 
lwin_fsd = 100;  threshold_fsd = 9;  shift_fsd = 100; lambda = 50;

% LTSD (window size, decision threshold and overlapping).
lwin_ltsd = 100;       threshold_ltsd = 9;   shift_ltsd = 10;

% 3) Get the decision signal of the FSD algorithm and the marker.
[V_fsd, T_fsd] = wag.fsd(input_signal, lwin_fsd, shift_fsd, 512, ...
    threshold_fsd);
[marker_fsd, T_fsd_expanded] = wag.compEstMark(V_fsd, T_fsd, input_signal, ...
    lwin_fsd, shift_fsd);

% 4) Get the decision signal of the LTSD algorithm and the marker.
[V_ltsd, T_ltsd] = wag.ltsd(input_signal, lwin_ltsd, shift_ltsd, 512, ...
    threshold_ltsd);
[marker_ltsd, T_ltsd_expanded] = wag.compEstMark(V_ltsd, T_ltsd, ...
    input_signal, lwin_ltsd, shift_ltsd);

if strcmpi(showPlotsCheck,'yes')    
figure
subplot(2, 1, 1)
plot(T_fsd_expanded)
hold on
plot(threshold_fsd * ones(1, length(T_fsd_expanded)), 'r')
legend('Detector output (FSD)', 'Detection threshold')
subplot(2, 1, 2)
plot(T_ltsd_expanded)
hold on
plot(threshold_ltsd * ones(1, length(T_ltsd_expanded)), 'r')
legend('Detector output (LTSD)', 'Detection threshold')

figure
subplot(2, 1, 1)
plot(input_signal)
hold on
plot(marker_fsd.*400, 'r')
legend('Input signal','FSD decision')
subplot(2, 1, 2)
plot(input_signal)
hold on
plot(marker_ltsd.*400, 'r')
legend('Input signal','LTSD decision')
end

% Determinate the initial and end point of each interval where we need to 
% find the peaks, i.e, the first activity period of each cycle. 
edges = find(diff(marker_ltsd)~=0);
initcross = edges(1:4:length(edges));
finalcross = edges(2:4:length(edges));

% Calculate the peaks in each interval.
for k = 1:length(initcross)                                 
 
    % Find all peaks greater than threshold in the interval specified by the 
    % minimum peak calculated above and the end of the interval.
    [peak_values_l, peak_locations_l] = findpeaks(g_Y_left_shank_1_C(...
                                        initcross(k):finalcross(k)));

    % Store the index of the highest positive peak .                                      
    sync_peaks_l(k) = find(...
                     g_Y_left_shank_1_C(initcross(k):finalcross(k)) ...
                     == max(peak_values_l), 1, 'first') + initcross(k) - 1;
 
    % Find all peaks greater than threshold in the interval specified by the 
    % minimum peak calculated above and the end of the interval.
    [peak_values_r, peak_locations_r] = findpeaks(g_Y_right_shank_1_C(...
                                        initcross(k):finalcross(k)));

    % Store the index of the highest positive peak .                                      
    sync_peaks_r(k) = find(...
                     g_Y_right_shank_1_C(initcross(k):finalcross(k)) ...
                     == max(peak_values_r), 1, 'first') + initcross(k) - 1;

end

% Evaluate if the patient steps with the left or right limb first for each
% cycle and store the first peak in sync_peaks, respectively, then sort it.
% In this case, we use the gyroscopo data.
sync_peaks_gyro = [sync_peaks_r(sync_peaks_l > sync_peaks_r), ...
              sync_peaks_l(sync_peaks_r > sync_peaks_l)];
sync_peaks_gyro = sort(sync_peaks_gyro);

% Detected sync-peaks.
if strcmpi(showPlotsGyroShank,'yes')
figure()
subplot(2, 1, 1);
plot(time_GW, g_Y_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l), g_Y_left_shank_1_C(sync_peaks_l), 'r.');

title('Angular Velocity  (Gyroscope) of the left shank with detected sync-peaks');
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');

subplot(2, 1, 2)
plot(time_GW, g_Y_right_shank_1_C, 'g');
hold on;
plot(time_GW(sync_peaks_r), g_Y_right_shank_1_C(sync_peaks_r), 'm.');

title('Angular Velocity (Gyroscope) of the right shank with detected sync-peaks');   
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');

% Comparation between peaks detection in Acc and Gyro signals.
figure ()
plot(time_GW(sync_peaks_acc), a_Z_right_shank_1_C(sync_peaks_acc),'m.');
hold on;
plot(time_GW(sync_peaks_gyro),g_Y_right_shank_1_C( sync_peaks_gyro), 'r.');
hx = graph2d.constantline(time_GW(sync_peaks_gyro), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'g');
changedependvar(hx,'x');
legend ('Acc', 'Gyro', 'Location', 'NorthEastOutside');
title('Comparation between peaks detection in Acc and Gyro signals');  

end
 
% Show other signals of the trunk for acceletometer data and the gyroscope
% data.
if strcmpi(showPlotsGyroTrunk,'yes')
figure()
subplot(3, 1, 1);
plot(time_GW, g_X_center_trunk_1_C);
hold on;
plot(time_GW(sync_peaks_gyro), g_X_center_trunk_1_C(sync_peaks_gyro), 'r.');

title('Orientation in the trunk (axe X) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Angle (deg)');

subplot(3, 1, 2)
plot(time_GW, g_Y_center_trunk_1_C);
hold on;
plot(time_GW(sync_peaks_gyro), g_Y_center_trunk_1_C(sync_peaks_gyro), 'r.');

title('Orientation in the trunk (axe Y) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Angle (deg)');

subplot(3, 1, 3)
plot(time_GW, g_Z_center_trunk_1_C);
hold on;
plot(time_GW(sync_peaks_gyro), g_Z_center_trunk_1_C(sync_peaks_gyro), 'r.');

title('Orientation in the trunk (axe Z) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Angle (deg)');

figure()
subplot(3, 1, 1);
plot(time_GW, a_X_center_trunk_3_C);
hold on;
plot(time_GW(sync_peaks_acc), a_X_center_trunk_3_C(sync_peaks_acc), 'r.');

title('Acceleration in the trunk (axe X) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Acc (g)');

subplot(3, 1, 2)
plot(time_GW, a_Y_center_trunk_3_C);
hold on;
plot(time_GW(sync_peaks_acc), a_Y_center_trunk_3_C(sync_peaks_acc), 'r.');

title('Acceleration in the trunk (axe Y) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Acc (g)');

subplot(3, 1, 3)
plot(time_GW, a_Z_center_trunk_3_C);
hold on;
plot(time_GW(sync_peaks_acc), a_Z_center_trunk_3_C(sync_peaks_acc), 'r.');

title('Acceleration in the trunk (axe Z) with detected sync-peaks ');
xlabel('Time in s');
ylabel('Acc (g)');

end
end

% Show a completion message to indicate a successful synchronisation of
% all patient data in the selected Excel-file.

% Read a picture to show in the box.
icon = imread('ok.jpg');

% Show completion message.
msgbox(['Synchronisation Completed! You can find the synchronised files in', ...
        '/APA/First iteration/data/Synchronised'], ...
        'Synchronisation Completed!', 'custom', icon);


