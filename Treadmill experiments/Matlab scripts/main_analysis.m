% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||| QUALISYS GAIT DATA ANALYSIS ||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% -------------------------------------------------------------------------

% The following script is used to analyze gait data gathered using a
% Qualisys optical system. 

% *************************************************************************
% - Authors: Dr. Prof. med. Kai B?tzel and Dr. Eng. Alberto Olivares.
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

% calculate angular position of each marker trialgle in 2D
% and put them into data array Q_leg_pitch
% [pitch_KF_right_shank, pitch_KF_right_thigh, pitch_KF_left_shank, pitch_KF_left_thigh];

Q_leg_pitch = zeros(length(data),4);


order = {   'right lo shank','right up shank','right back shank';...
            'right lo thigh','right up thigh','right back thigh';...
            'left lo shank','left up shank','left back shank';...
            'left lo thigh','left up thigh','left back thigh';...
            'right hip','back hip','left hip'};
x=1;z=3;lag = zeros(4,1);
figure(1)
%  big_ax=axes('position', [0 0 1 1],'visible', 'off');

for n=1:4                                                                       % 4 leg segments
    m1=0;m2=0;m3=0;
    
    for n2=1:3                                                                  % 3 markers
        ind = 1; %find(strcmp(order{n,n2},data_labels));
        if isempty(ind)
            fprintf('Error: cannot find match for %s\n',order{n,n2})
            beep
        else
            switch n2
                case 1, m1 = ind;
                case 2, m2 = ind;
                case 3, m3 = ind;
            end
        end
    end
    
    
    
    seg_1 = -atan(squeeze((data(m2,x,:)-data(m1,x,:))./(data(m1,z,:)-data(m2,z,:))));        % marker 1,2: x/z
    seg_2 = -atan(squeeze((data(m2,z,:)-data(m3,z,:))./(data(m2,x,:)-data(m3,x,:))));        % marker 2,3: z/x
    seg_3 = -atan(squeeze((data(m1,z,:)-data(m3,z,:))./(data(m1,x,:)-data(m3,x,:))));        % marker 3,1: x/z

     Q_leg_pitch(:,n) = mean([seg_1 seg_2 seg_3],2);
     
end


