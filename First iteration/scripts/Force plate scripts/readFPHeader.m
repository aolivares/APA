function [header, first_line, samp_time, count, col_first_block] = ...
    readFPHeader(filename)

% FUNCTION readFPHeader Reads the header of the data files gathered by
% the Zebris FDM-S Multifunction Force-measuring Plate.
%
% - Input:
%    |_ 'filename': Path to the data file (string).
%
% - Output:
%    |_ 'header': Structure with the following components:
%       |_ 'header{1}': Name of the patient or patient identifier (string).
%       |_ 'header{2}': Name of the application, i.e. type of experimet and
%                       type of data file (string).
%       |_ 'header{3}': Data labels i.e. name of the data columns (struct).
%       |_ 'header{4}': Calibration flag (int).
%       |_ 'header{5}': Creation date. Format: dd/mm/yy hh:mm:ss (string).
%   |_ 'first_line':      Number of the row in which the data start (int).
%   |_ 'samp_time':       Sampling time in miliseconds (int).
%   |_ 'count':           Number of data blocks (int).
%   |_ 'col_first_block': Number of columns of the first block (int).
%
% -------------------------------------------------------------------------
% * Authors:      - Prof. Dr. Med. Kai B�tzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Ver�nica Torres (2): 
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
% * Last modification: 08/10/2014
% -------------------------------------------------------------------------

% Open data file.
fid = fopen(filename);

% Show name of data file.
disp(filename);

% Read the first 19 lines (this is where the header info is).
C = textscan(fid, '%s', 19, 'Delimiter', '\n', 'whitespace', '');   

% Read name (or patient identifier) of the patient.
ind = find(strncmp('Per', C{1}, 3));                          
if ~isempty(ind)
        nn = textscan(char(C{1}(ind)), '%s');
        nn = char(nn{:}); 
        nn(1, :) = []; 
        [r, c] = size(nn);                    
        header{1} = reshape(nn', 1, r*c);
end

% Read the name of the applicataion (type of experiment and data file).
ind = find(strncmp('App', C{1}, 3));                          
nn = textscan(char(C{1}(ind)), '%s%s');
header{2} = char(nn{2});

% Read the file creation date.
ind = find(strncmp('Cre', C{1}, 3));                          
nn = textscan(char(C{1}(ind)), '%s%s%s%s');
header{5} = [char(nn{3}), ' ', char(nn{4})];

% Read sampling time in ms.
ind = find(strncmp('Fre', C{1}, 3));                          
nn = textscan(char(C{1}(ind)), '%s%n');
samp_time = 1000 / nn{2}; 

% Read number of data blocks.
ind = find(strncmp('Cou', C{1}, 3));                           
nn = textscan(char(C{1}(ind)), '%s%n');
count = nn{2};

% Read data labels and number of columns in the data blocks.
ind = find(strncmp('Tim', C{1}, 3));                           
if strcmp(header{2}, 'trplatf') || strcmp(header{2}, 'statplat')
     data_labels = textscan(char(C{1}(ind)), '%s', 'delimiter', '\t');
     col_first_block = length(data_labels{1});                                  
     header(3) = data_labels;                                             
end

% Read calibration flag.
ind = find(strncmp('Cal', C{1}, 3));                           
if strcmp(header{2}, 'trplatf')
    % Get calibration data of the Zebris Armswing marker.
    calib = textscan(char(C{1}(ind)), '%n', 'delimiter', '\t', ...
        'treatAsEmpty', 'Calibr');  
    header{4} = calib{1}(~isnan(calib{1}));
else
    header{4} = 0;
end

% Get the number of the row where the data start.
first_line = ind + 1; 

% Close data file.
fclose(fid);

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF READ_FP_HEADER FUNCTION
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\