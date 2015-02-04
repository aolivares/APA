% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------------ Extract Signals --------------------------------
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
% * Last modification: 04/02/2015
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
% * 1) Selection of the excel file and obtain the GW file names.
% 
% * 2) Extraction FP data from  .txt files.
% 
% * 3) Calibrate GW data.
% 
% * 4) Synchronisation.
% 
% * 5) Save the synchronised data in *.mat file for each person.
% -------------------------------------------------------------------------

%% 

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

% -------------------------------------------------------------------------
% 1) Select, read and obtain information from the excel file.
% -------------------------------------------------------------------------

% Select only one data files with a dialog box (only .xlsl files).
[filename_excel, filepath] = uigetfile('*.xlsx', ...
    'Select the Excel file with tha data (.xlsx)', '../../data');

% Read data from *.xlsx where are stored all filenames and other 
% interesting information.
[~,file_excel] = xlsread(['../../data/' filename_excel]);
[rows,columns] = size(file_excel);

% The two first rows are unuseful so we obtain the data beginning from 
% second rows. We extract GW names that are in the third column, the 
% number of FP files and their names to read them afterwards.
ind = 1;

for i=15:rows-1
        
        % Check if the cell is empty
        if  isempty( char(strtrim(file_excel(i,4))) ) == 0  

            % We extract the GW name that identifies the patient.
            filename_GW_total (ind) = strtrim(file_excel(i,4)); 
            index_GW (ind)=i;
            ind= ind+1;
        end
        
end

index_GW= [index_GW rows+1];

%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Data processing for each patient identified previously.
% -------------------------------------------------------------------------
%--------------------------------------------------------------------------

for j=1:length(index_GW)-1
    
% -------------------------------------------------------------------------
% 2) Extraction FP data from  .txt files.
% -------------------------------------------------------------------------
    start_data_patient=index_GW(j);
    end_data_patient=index_GW(j+1)-1;
    
    filename_complete_FP=strtrim(file_excel(...
        start_data_patient:end_data_patient,1));
    
    n_files = length(filename_complete_FP);
    
    % ---------------------------------------------------------------------
    % 2.1) Read the blocks of data and compute the interesting signals.
    % ---------------------------------------------------------------------

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
            readFPHeader(fullfile('../../data/ForcePlate/Raw/',filename_FP));

        % Read formatted data from text file. As always "fid" is a file 
        % identifier that we obtain with "fopen".
        fid = fopen(fullfile('../../data/ForcePlate/Raw/', filename_FP));

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

        % Values: (1): COP (Center of Gravity) X (right-left).
        %         (2): COP Y (anteroposterior i.e. front-to-back).
        %         (3): Posterior margin (only Y).
        %         (4): Anterior margin (only Y is loaded).
        %         (5): Overall pressure.

        % Define the variable that contains the force data.
        force_cell = zeros(lines, columns);
        force_cell_complete_cycle = zeros(n, lines, columns);

        fid = fopen(fullfile('../../data/ForcePlate/Raw/',filename_FP));
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

            % Index for data
            c = b - data_start + 1; 

            % Read 8 lines which contain some irrelevant information
            CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
            C3 = textscan(fid, format, lines);
            for cc = 1:columns
               force_cell(:, cc) = C3{cc + 1};
            end

            % Store the matriz with the force of each cell in the correct
            % position.
            force_cell_complete_cycle(n,:,:) = force_cell;

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
                m = sum(force_cell(:, an:en));   

                % Check whether any data is >0.
                if any(m)    
                    % Put the COP X in the rigth position. We multiply this  
                    % result by 8.5 because each cell has 8.5 mm (as much width
                    % as heigth)
                    lateral_COP(c,side) = 8.5.*sum(m .* (an:en)) ./ sum(m);

                    % Calculate COP Y.
                    m = sum(force_cell(:, an:en), 2);                                          
                    antpost_COP_cycle(c, side) = 8.5.*sum(m .* (1:lines)') ...
                        ./sum(m);

                    % Calculate rearmost foot pressure point (heel position) 
                    % (y only)
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
        time_FP{i} = time_FP_cycle;
        force_cells{i} = force_cell_complete_cycle;
        ML_COP{i} = medlateral_COP_cycle';
        AP_COP{i} = antpost_COP_cycle';
        force_sum{i} = force_complete_cycle';
        midline{i} = midline_mm;

    end

   clearvars -except force_sensors time_FP force_cells ML_COP ...
    AP_COP force_sum midline filename_GW_total index_GW file_excel j

    
    
    
    
    
% -------------------------------------------------------------------------
% 3) Calibrate GW data.
% -------------------------------------------------------------------------  
   
    filename_GW= filename_GW_total(j);
    
  
    
 
% -------------------------------------------------------------------------
% 4) Synchronisation.
% -------------------------------------------------------------------------
    

% -------------------------------------------------------------------------
% 5) Save the synchronised data in *.mat file for each person.
% -------------------------------------------------------------------------
    
    
    
    
end



% -------------------------------------------------------------------------
% 6) Show a completion message to warn the end of the run.
% -------------------------------------------------------------------------

% Read a picture to show in the box.
icon = imread('ok.jpg');

% Show the massage to warn the end.
msgbox('Operation Completed!!!!', 'Success', 'custom', icon);


