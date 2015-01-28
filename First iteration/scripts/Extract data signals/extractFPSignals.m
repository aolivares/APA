
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------- Extract FP Signals --------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai B?tzel (1): 
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
% * Last modification: 26/01/2015
% -------------------------------------------------------------------------
% INFORMATION: This file contains the routine to extract the FP signals 
% from the ten separates forceplate cycles. After, these signals are stored 
% in a FP_data_*.mat file.
% The stored output file includes seven  cell type variables:
% 
% * force_sensors       : four signals that contains pressure of the four
%                          forceplate sensors.
% * time_FP             : time signal.
% * force_cell_complete : matrices collection with the pressure of every
%                         forceplate cell computed for each time frame.
% * medlateral_COP      : medio-lateral COP (Center of Pressure) for three
%                         differents cases organized by rows: rigth foot, 
%                         left foot and both feet  respectively.
% * antpost_COP         : antero-posterior COP (Center of Pressure) for 
%                         three differents cases organized by rows: rigth  
%                         foot, left foot and both feet  respectively.
% * force_complete      : overall force for  three differents cases   
%                         organized by rows: rigth foot, left foot and both 
%                         feet  respectively.
% * midline             : midline between both feets.
% 
% Each of them contains the signals of all cycles, organized by rows.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0) Clear workspace and close all figures.
% -------------------------------------------------------------------------
clear all; close all; clc;

% -------------------------------------------------------------------------
% 1) Obtain the files names.
% -------------------------------------------------------------------------

% Select data files with a dialog box (only .txt files). Press 'Ctrl' and
% click all data of the same patient.
[filename_complete, filepath] = uigetfile('*.txt', ...
    'Select Force Plate data file (.txt)', '../../data/ForcePlate/Raw',...
    'MultiSelect','on');

% Determine the number of files that we have selected.
[aux,n_files] = size(filename_complete);

% Check if filename_complete is a cell because if it isn't like that, only 
% one file has been selecting.
tf = iscell(filename_complete);
if tf == 0
    msgbox('Please, you must select at least two files  ',...
    'Error','error');
    beep;
    
end

% -------------------------------------------------------------------------
% 3) Sort the selected names.
% -------------------------------------------------------------------------

% Sort the selected names for they are in the correct order to match with
% the appropiate GW signals.
[filename_complete,I] = sort(filename_complete);

% Rearrange the names again, because somo files like ...03 2 must be after 
% ..03. We differenciate this name file name.
ind = 0;

for i = 1:n_files
    
    % Compare the length of one namefile and the next. If the first is
    % shorter than the second, we change tha values to put the names in the
    % correct order.
    ind = ind+1;
    
    % Check we don't access non-existen array position.
    if(ind < n_files)
        if length(filename_complete{ind}) > length(filename_complete{ind+1})
            
          % We use a auxiliar variable to not overwrite the value of the 
          % interchanged namefile
          aux_filename = filename_complete{ind};
          filename_complete{ind} = filename_complete{ind+1};
          filename_complete{ind+1} = aux_filename;
          ind = ind+1;
          
        end
    end
end

% -------------------------------------------------------------------------
% 4) Read the blocks of data and compute the interesting signals.
% -------------------------------------------------------------------------

% Define the variables used to store the data.
force_sensors = cell(n_files, 1);
time_FP = cell(n_files, 1);
force_cell_complete = cell(n_files, 1);
medlateral_COP = cell(n_files, 1);
antpost_COP = cell(n_files, 1);
force_complete = cell(n_files, 1);
midline = cell(n_files, 1);

% Read each selected data file.
for i = 1:n_files
    
    % Convert the namefile in a char type to be able to read the file.
    filename = char(filename_complete(i));

    % Read data file.
    [header, first_line, abta, count, col_first_block] = ...
        readFPHeader(fullfile(filepath,filename));

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
    force_cell_complete{i} = force_cell_complete_cycle;
    medlateral_COP{i} = medlateral_COP_cycle';
    antpost_COP{i} = antpost_COP_cycle';
    force_complete{i} = force_complete_cycle';
    midline{i} = midline_mm;
    
end

% -------------------------------------------------------------------------
% 5) Save all data.
%--------------------------------------------------------------------------

% Obtain the name for the data that corresponding to the patient.
name_file = textscan(filename,'%s','Delimiter',',');
name_file = char(name_file{1}{1});
name_file = strcat('FP_data_',name_file);

% Save all data.
save(['../../data/ForcePlate/Preprocessed/' name_file '.mat'],  ...
    'force_sensors','time_FP', 'force_cell_complete', 'medlateral_COP', ...
    'antpost_COP', 'force_complete','midline');

% -------------------------------------------------------------------------
% 6) Show a completion message to warn the end of the run.
% -------------------------------------------------------------------------

% Read a picture to show in the box.
icon = imread('ok.jpg');

% Show the massage to warn the end.
msgbox('Operation Completed!!!! You can find the created file in ForcePlate folder!',...
    'Success', 'custom', icon);

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF EXTRACT FP SIGNALS FILE
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\