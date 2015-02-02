
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
% * Last modification: 29/01/2015
% -------------------------------------------------------------------------
% INFORMATION: This file contains the routine to extract the FP signals 
% from the ten separates forceplate cycles. After, these signals are stored 
% in a FP_data_*.mat file.
% The stored output file includes seven  cell type variables:
% 
% * force_sensors       : four signals that contains pressure of the four
%                          forceplate sensors.
% * time_FP             : time signal.
% * force_cells         : matrices collection with the pressure of every
%                         forceplate cell computed for each time frame.
% * ML_COP              : medio-lateral COP (Center of Pressure) for three
%                         differents cases organized by rows: rigth foot, 
%                         left foot and both feet  respectively.
% * AP_COP              : antero-posterior COP (Center of Pressure) for 
%                         three differents cases organized by rows: rigth  
%                         foot, left foot and both feet  respectively.
% * force_sum           : overall force for  three differents cases   
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

% Select only one data files with a dialog box (only .txt files).
[filename_one_file, filepath] = uigetfile('*.txt', ...
    'Select Force Plate data file (.txt)', '../../data/ForcePlate/Raw');

% Obtain the data name that corresponding to the patient to compare
% with the rest of files names stored in .xlsx.
name_file = textscan(filename_one_file,'%s','Delimiter',',');
name_file = char(name_file{1}{1});

% -------------------------------------------------------------------------
% 2) Obtain all name files belonging to the same person.
% -------------------------------------------------------------------------

% Read data from APA_Munic filenames.xlsx where are stored all file
% names and other interesting information.
sheet = 'APA_liste_3';
[~,txt] = xlsread('../../data/APA_Munic filenames.xlsx',sheet);
[rows,columns]=size(txt);

% The four first rows are unuseful. We obtain only the FP names that are in
% the first column.
ind=1;

for i=5:rows-1
    
    % We extract the name that identifies the patient to compare with the 
    % name file selected previously.
    name_individual_file = strtrim(txt(i,1));               
    name_individual_file = textscan(name_individual_file{1,1},'%s',...
        'Delimiter',',');
    name_individual_file = char(name_individual_file{1}{1});
    
    % Check if the string lengths are the same because you can't compare 
    % char arrays with different sizes.
    if length( name_individual_file) == length(name_file )
        
        % Check if the file names are the same, i.e if the file name 
        % corresponds with person selected. These names are stored in a
        % variable to access data file afterwards.
        if (name_individual_file == name_file)
            filename_complete( ind) = strtrim(txt(i,1));
            ind = ind+1;
        end
    end

end

n_files = length(filename_complete);

% -------------------------------------------------------------------------
% 3) Read the blocks of data and compute the interesting signals.
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
    force_cells{i} = force_cell_complete_cycle;
    ML_COP{i} = medlateral_COP_cycle';
    AP_COP{i} = antpost_COP_cycle';
    force_sum{i} = force_complete_cycle';
    midline{i} = midline_mm;
    
end

% -------------------------------------------------------------------------
% 4) Save all data.
%--------------------------------------------------------------------------

% Create the output file name.
 name_file = strcat('FP_data_',name_file);

% Save all data.
save(['../../data/ForcePlate/Preprocessed/' name_file '.mat'],  ...
    'force_sensors','time_FP', 'force_cells', 'ML_COP', ...
    'AP_COP', 'force_sum','midline');

% -------------------------------------------------------------------------
% 5) Show a completion message to warn the end of the run.
% -------------------------------------------------------------------------

% Read a picture to show in the box.
icon = imread('ok.jpg');

% Show the massage to warn the end.
msgbox('Operation Completed!!!! You can find the created file in ForcePlate folder!',...
    'Success', 'custom', icon);

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF EXTRACT FP SIGNALS FILE
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\