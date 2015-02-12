% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||| QUALISYS GAIT DATA ANALYSIS ||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% -------------------------------------------------------------------------

% The following script is used to analyze gait data gathered using a
% Qualisys optical system. 

% *************************************************************************
% - Authors: Dr. Prof. med. Kai Bötzel and Dr. Eng. Alberto Olivares.
% - Entities: Ludwig-Maximilian University and University of Granada.
% - Last revision: 10/03/2014.
% *************************************************************************

% -------------------------------------------------------------------------
% 0) Initial configuration.
% -------------------------------------------------------------------------
clear all;

% -------------------------------------------------------------------------
% 1) Load data.
% -------------------------------------------------------------------------

% Prompt to select the data file which is to be analyzed. 
[filename, filepath] = uigetfile('*.mat', 'Select the QS data file');

% Load data.
load(fullfile(filepath, filename));

% Remove the '.mat' extension from the file name.
filename = filename(1:end-4);

% -------------------------------------------------------------------------
% 2) Read data.
% -------------------------------------------------------------------------

% Read sampling rate. The sampling rate is stored inside the 'FrameRate'
% field of structure which name is equal to the value of 'filename'.
SampRate = eval([filename, '.FrameRate']);

% Read markers data. Data are stored in the 'Trajectories.Labeled.Data'
% field
data = eval([filename, '.Trajectories.Labeled.Data']);

% Remove the marker error information which is also stored by the Qualisys
% software. By doing this we only leave the (X, Y, Z) data.
data(:,4,:) = []; 

% Change the polarity of the data to match the desired coordinate frame.
data = data * -1;    

% Get size of the data matrix.
[n_markers, dummy, n_frames] = size(data);

% Extract data labels.
data_labels = eval([filename,'.Trajectories.Labeled.Labels']);

% Apply some slight modifications to the name of the labels. More
% especifically the underscores are removed and the substrings 'lo' and
% 'up' are replaced by 'lower' and 'upper' respectively.
for n = 1:n_markers
    data_labels{n} = strrep(data_labels{n},'_',' ');
    if ~isempty(strfind(data_labels{n},'lo'))
       data_labels{n} = strrep(data_labels{n},'lo','lower'); 
    end
    if ~isempty(strfind(data_labels{n},'up'))
       data_labels{n} = strrep(data_labels{n},'up','upper'); 
    end
end     

% Remove data structure for memory reasons.
eval(['clear ' filename]); 

% -------------------------------------------------------------------------
% 3) Interpolate data to remove NaN values (if existing).
% -------------------------------------------------------------------------
% The data matrix may have some 'NaN' values which are written by the
% Qualisys software in those instants in which the markers are not visible
% by any of the eight cameras. Markers are mainly ocluded by loose clothes
% so it is necessary to deal with this situation. We will apply
% interpolation to solve this issue.

% We first show how many data are missing for each one of the markers.
nan_list = cell(1,n_markers);
for n = 1:n_markers
    nan_list{n} = [num2str(n),')  ', 'Marker: ',data_labels{n},'. Number of Nan values:', ...
        num2str(sum(isnan(squeeze(data(n,1,:)))))];
end

[Selection,ok] = listdlg('ListString',nan_list,'Name',...
    'Select the marker data you wish to correct','ListSize',[450 400],'SelectionMode',...
    'multiple');


