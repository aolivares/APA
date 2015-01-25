
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------- Extract FP Signals --------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Bötzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Verónica  Torres (2): 
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
% * Last modification: 20/01/2015
% -------------------------------------------------------------------------
% INFORMATION: This file contain the routine to extract the FP signals from 
% the ten separates forceplate cycles. After, these signals are stored in a 
% FP_data_*.mat file.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clear workspace and close all figures.
% -------------------------------------------------------------------------
clear all; close all; clc;

% -------------------------------------------------------------------------
% Read the blocks of data
% -------------------------------------------------------------------------

% Select data files with a dialog box (only .txt files). Press 'Ctrl' and
% click all data of the same patient.
[filename_complete, filepath] = uigetfile('*.txt', ...
    'Select all forceplate data file of the same patient(.txt)', ...
    '../../data/ForcePlate','MultiSelect','on');

% Check if filename_complete is a cell because if it isn't like that, only 
% one file has been selecting.
tf=iscell(filename_complete);
if tf==0
    msgbox('Please, you must select at least two files  ',...
    'Error','error');
    beep;
    
end


% -------------------------------------------------------------------------
% Rearrange all selected files to match with GW signals correctly.
% -------------------------------------------------------------------------

% Determine the number of files that we have selected.
[x, n_files]=size(filename_complete);
[filename_complete,I]=sort(filename_complete);

% -------------------------------------------------------------------------
% Data processing.
% -------------------------------------------------------------------------

% Use the first name to create el completed name of the output file.
filename=char(filename_complete(1));

% Obtain the name for the data that corresponding to the patient.
name_file=textscan(filename,'%s','Delimiter',',');
name_file=char(name_file{1}{1});
name_file=strcat('FP_data_',name_file);

% Read each selected data file.
for i=1:n_files
    
    % Convert the namefile in a char type to be able to read the file.
    filename=char(filename_complete(i));

    % Read data file.
    [header, first_line, abta, count, col_first_block] = ...
        readFPHeader(fullfile(filepath,filename));
    
    % Read first block of data. The first block contains 5 columns:
        %   (1) Time (ms).
        %   (2) 'Li Vorfuß,N':  Front left foot pressure sensor.
        %   (3) 'Li Rückfuß,N': Back left foot pressure sensor.
        %   (4) 'Re Vorfuß,N':  Front right foot pressure sensor.
        %   (5) 'Re Rückfuß,N': Back right foot pressure sensor.

    % Read formatted data from text file. As always "fid" is a file 
    % identifier that we obtain with "fopen".
    fid = fopen(fullfile(filepath, filename));

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
    
    % Define the four pressure signals of the diferents platforce sensors 
    % as well as the time.
    time_FP=data(data_start:data_end, 1);
    left_front_foot_force=data(data_start:data_end, 2);
    left_back_foot_force=data(data_start:data_end, 3);
    right_front_foot_force=data(data_start:data_end, 4);
    right_back_foot_force=data(data_start:data_end, 5);


    % Read 10 lines which contain information concerning the following data 
    % blocks.
    C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');   

    % Find out no of columns per block.
    nn = textscan(char(C2{1}(5)), '%s%n');                              
    columns = nn{2};

    % Find out no of lines per block 
    nn = textscan(char(C2{1}(6)), '%s%n');                             
    lines = nn{2};

    % Set format to read the data, first value is a 'string' variable, the 
    % rest are 'double' variables.
    format = ['%s', repmat('%n', 1, columns)];                            

    % Read data block by block and evaluate only those which contain data
    % only to define the midline. 
    m = zeros(columns,1);
    
    % 

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
    % data_file = zeros(n,3,4); 

    % Values: (1): COP (Center of Gravity) X (right-left).
    %         (2): COP Y (anteroposterior i.e. front-to-back).
    %         (3): Posterior margin (only Y).
    %         (4): Anterior margin (only Y is loaded).
    %         (5): Overall pressure.

    % Define the variable that contains the force data of each cell.
    force_cell = zeros(lines, columns);
    
    fid = fopen(fullfile(filepath,filename));
    dummy_data = textscan(fid, '%n', col_first_block * count, 'headerlines',...
        first_line - 1);
    C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');

    % Skip these data they contain no values
    for b = 1:data_start - 1  
         % Read 8 lines which contain some irrelevant information
         CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
         C3 = textscan(fid, format, lines);
    end


    for b = data_start:data_end
        % Index for data3
        c = b - data_start + 1; 
        % Read 8 lines which contain some irrelevant information
        CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
        C3 = textscan(fid, format, lines);
        for cc = 1:columns
           force_cell(:, cc) = C3{cc + 1};
        end

        % Calculate COP in the differents cases.
        % Side 1 right, side 2 left, side 3 both
        for side = 1:3                                                              
            switch side
                case 1   
                 % Column vector. Calculate COP X to represent the markers.
                     m = sum(force_cell(:, 1:ml));   

                 % Check whether any data is >0.
                    if any(m)    
                        % Put the COP X in the rigth position. We multiply 
                        % this result by 8.5 because each cell has 8.5 mm  
                        % (as much width as heigth)
                        right_lateral_COP(c) = 8.5.*sum(m .* (1:ml)) ./ ...
                            sum(m);
                        
                        % Definition of the midline in mm.
                        ml_mm=8.5*ml;
                        
                        % Calculate COP X with respect to midline.
                        right_lateral_COP=right_lateral_COP - ml_mm;
                        
                        % Define the 0 number like not a number.
                        right_lateral_COP(right_lateral_COP==0) = NaN;

                        % Calculate COP Y.
                        m = sum(force_cell(:, 1:ml), 2);                                          
                        right_antpos_COP (c) = 8.5.*sum(m .* (1:lines)') ./ ...
                            sum(m);
                        
                        % Define the 0 number like not a number.
                        right_antpos_COP(right_antpos_COP==0) = NaN;
                        
                        % Calculate rearmost foot pressure point 
                        % (heel position,y only)
                        right_force (c) = sum(m);

                        % Define the 0 number like not a number.
                        right_force(right_force==0) = NaN;
                    end
                    
                case 2  
                 % Column vector. Calculate COP X to represent the markers.
                     m = sum(force_cell(:, ml+1:columns));   

                 % Check whether any data is >0.
                    if any(m)    
                        % Put the COP X in the rigth position. We multiply 
                        % this result by 8.5 because each cell has 8.5 mm  
                        % (as much width as heigth)
                        left_lateral_COP(c)=8.5.*sum(m .* (ml+1:columns))...
                            ./ sum(m); 
                        
                        % Define the 0 number like not a number.
                        left_lateral_COP(left_lateral_COP==0) = NaN;
                        
                        % Definition of the midline in mm.
                        ml_mm=8.5*ml;
                        
                        % Calculate COP X with respect to midline.
                        left_lateral_COP=left_lateral_COP - ml_mm;

                        % Calculate COP Y.
                        m = sum(force_cell(:, ml+1:columns), 2);                                          
                        left_antpos_COP (c) = 8.5.*sum(m .* (1:lines)') ./ ...
                            sum(m);
                        
                        % Define the 0 number like not a number.
                        left_antpos_COP(left_antpos_COP==0) = NaN;

                        % Calculate rearmost foot pressure point 
                        % (heel position,y only)
                        left_force(c) = sum(m);
                        
                        % Define the 0 number like not a number.
                        left_force(left_force==0) = NaN;

                    end
                    
                case 3
                 % Column vector. Calculate COP X to represent the markers.
                     m = sum(force_cell(:, 1:columns));   

                 % Check whether any data is >0.
                    if any(m)    
                        % Put the COP X in the rigth position. We multiply 
                        % this result by 8.5 because each cell has 8.5 mm  
                        % (as much width as heigth)
                        both_lateral_COP (c)= 8.5.*sum(m .* (1:columns))...
                            ./ sum(m); 
                        
                        % Define the 0 number like not a number.
                        both_lateral_COP(both_lateral_COP==0) = NaN;
                        
                        % Definition of the midline in mm.
                        ml_mm=8.5*ml;
                        
                        % Calculate COP X with respect to midline.
                        both_lateral_COP=both_lateral_COP - ml_mm;

                        % Calculate COP Y.
                        m = sum(force_cell(:, 1:columns), 2);                                          
                        both_antpos_COP (c) = 8.5.*sum(m .* (1:lines)')...  
                            ./ sum(m);
                        
                        % Define the 0 number like not a number.
                        both_antpos_COP(both_antpos_COP==0) = NaN;

                        % Calculate rearmost foot pressure point 
                        % (heel position,y only)
                        both_force (c) = sum(m);
                        
                        % Define the 0 number like not a number.
                        both_force(both_force==0) = NaN;

                    end
            end
            
        end

    end

    % Close data file.
    fclose(fid);
    
    
    % Store in the data struct.
    FP_file=struct( 'time_FP', time_FP, 'left_front_foot_force',...
      left_front_foot_force, 'left_back_foot_force', left_back_foot_force,...
      'right_front_foot_force', right_front_foot_force,...
      'right_back_foot_force', left_front_foot_force', 'force_cell', ...
      force_cell, 'right_lateral_COP', right_lateral_COP, ...
      'right_antpos_COP',right_antpos_COP, 'right_force', right_force,...
      'left_lateral_COP', left_lateral_COP, 'left_antpos_COP', ...
      left_antpos_COP, 'left_force', left_force, 'both_lateral_COP', ...
      both_lateral_COP, 'both_antpos_COP', ...
      both_antpos_COP, 'both_force', both_force, 'midline', ml_mm, ...
      'data_start', data_start, 'data_end', data_end );
  
    % Fix a new name for all signal of one file
    FP_file_name=strcat('FP_data', num2str(i));
    eval([ FP_file_name ' = FP_file;' ]);
     
    % Save all data.
    if i==1
         save(['../../data/ForcePlate/' name_file '.mat'], FP_file_name);
    else
        save(['../../data/ForcePlate/' name_file '.mat'], FP_file_name,...
            '-append');
    end
end

% -------------------------------------------------------------------------
% Show a completion message to warn the end of the run.
% -------------------------------------------------------------------------

% fprintf('\nThe run is finished!!!!!!\n\nIt has created %s.mat\n',name_file);
% disp('You can see in ForcePlate folder')

icon=imread('ok.jpg');
msgbox('Operation Completed!!!! You can find the created file in ForcePlate folder!',...
    'Success','custom',icon);

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF EXTRACT FP SIGNALS FILE
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        