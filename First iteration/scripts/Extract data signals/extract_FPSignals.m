
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
[filename_Complete, filepath] = uigetfile('*.txt', ...
    'Select Force Plate data file (.txt)', '../../data/ForcePlate',...
    'MultiSelect','on');

% Determine the number of files that we have selected.
[aux,n_files]=size(filename_Complete);

% Create the variable that will contain all the results.
% FP_data=cell(n_files);

% Read each selected data file.
for i=1:n_files
    
    % Convert the namefile in a char type to be able to read the file.
    filename=char(filename_Complete(i));

    % Read data file.
    [header, first_line, abta, count, col_first_block] = ...
        read_fp_header(fullfile(filepath,filename));

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
    data_file = zeros(n,3,4); 

    % Values: (1): COP (Center of Gravity) X (right-left).
    %         (2): COP Y (anteroposterior i.e. front-to-back).
    %         (3): Posterior margin (only Y).
    %         (4): Anterior margin (only Y is loaded).
    %         (5): Overall pressure.

    % Define the variable that contains the force data.
    data_force = zeros(lines, columns);

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
           data_force(:, cc) = C3{cc + 1};
        end

        % Calculate COP in the differents cases.
        % Side 1 right, side 2 left, side 3 both
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
            m = sum(data_force(:, an:en));   

            % Check whether any data is >0.
            if any(m)    
                % Put the COP X in the rigth position. We multiply this 
                % result by 8.5 because each cell has 8.5 mm (as much width
                % as heigth)
                data_file(c,1,side) = 8.5.*sum(m .* (an:en)) ./ sum(m); 


                % Calculate COP Y.
                m = sum(data_force(:, an:en), 2);                                          
                data_file(c, 2, side) = 8.5.*sum(m .* (1:lines)') ./ sum(m);

                % Calculate rearmost foot pressure point (heel position) 
                % (y only)
                data_file(c, 3, side) = sum(m);

            end
          
        end

    end
    
    % Complete each block of data with the appropiate time. It's the same 
    % in all cases.
    data_file(:,1,4)=data(data_start:data_end, 1);
    data_file(:,2,4)=data(data_start:data_end, 1);
    data_file(:,3,4)=data(data_start:data_end, 1);

    % Close data file.
    fclose(fid);
    
    % Define the 0 number like not a number.
    data_file(data_file==0) = NaN;

    % Definition of the midline in mm
    ml_mm=8.5*ml;

    % Calculate COP X with respect to midline.
    data_file(:, 1, :)=data_file(:, 1, :)-ml_mm;
    
    % Store in the data cell.
    FP_data{i}=data_file;
    
end

% Obtain the name for the data that corresponding to the patient.
name_file=textscan(filename,'%s','Delimiter',',');
name_file=char(name_file{1}{1});
name_file=strcat('FP_data_',name_file);

% Save all data.
save(['forcePlate_data/' name_file '.mat'], 'FP_data');


% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF EXTRACT FP SIGNALS FILE
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        