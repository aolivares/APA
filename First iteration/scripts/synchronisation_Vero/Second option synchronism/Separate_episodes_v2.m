function [ signals_separates, peak_loc_GW_Pral ] = Separate_episodes( GW_signal)

% FUNCTION SEPARATE_EPISODES Divide a signal in several smaller signals
% that you choose manually and interpole them, synchronized previously
% using like reference the first peak off of the signal.
%
% Input parameters:
% |_ 'GW_signal':  full signal from GW, it contain all repeats.
%
% Output parameters:
% |_ 'signals_separates': matrix where each row is a signal asociated each
%                         repeat. 
% |_ 'peak_loc_GW_Pral' : location of peak off used to synchronize.
%

% *************************************************************************
% - Author: Alberto Olivares, Kai Bötzel, Rob Weiss and Veronica Torres.
% - Entity: University of Granada, Spain. 
% - Last revision: 08/21/2013.
% *************************************************************************

% -------------------------------------------------------------------------
% Import GaitWatch functions library.
% -------------------------------------------------------------------------

gw = gwLibrary;

% -------------------------------------------------------------------------
% Separate the differents episodes 
% -------------------------------------------------------------------------

% Select the initial and final point of each repeat manually.
index= gw.getDCindexes (GW_signal,...
    'SELECT THE INITIAL AND FINAL POINTS OF EACH EPISODES');
close all;

% Set the number of repeats.
numb_episodes=length(index)-1;

% Determine the logest interval.
% Build a vector that contain the lenght of all intervals.
for i=1:numb_episodes-1
    gap (i)=index(i+1)-index(i);
end

% Select the longest interval that we`ll use like reference to interpolate 
% the rest of signals.
[maximum,index_max]=max(gap);
signalPral=GW_signal(index(index_max):index(index_max+1));


% Find the off-peak in GW  pral signal.
% Invert the signal to be able to use findpeaks cause we need to detect
% a minimum in the signal.

[peak_value, peak_location] = findpeaks(signalPral,'minpeakheight',1.4);
peak_loc_GW_Pral=peak_location(1);

% -------------------------------------------------------------------------
% Interpolte and obtain the rest of signals.
% -------------------------------------------------------------------------

% Initialize the output matrix that will contain all separated signals.
signals_separates=zeros(numb_episodes,maximum+1);


for i=1:numb_episodes
    
    % Obtain each interval.
    signal=GW_signal(index(i):index(i+1));
    
    % Find the off-peak in GW signal

    [peak_value, peak_location] = findpeaks(signal,'minpeakheight',1.4);
    peak_loc_GW=peak_location(1);

    % Interpolation
    div=peak_loc_GW/peak_loc_GW_Pral;
    xi1=1:div:peak_loc_GW;
    signal1=interp1(1:peak_loc_GW,signal(1:peak_loc_GW),...
        xi1,'spline');

    div=(length(signal)-peak_loc_GW)/(length(signalPral)-peak_loc_GW_Pral);
    xi2=peak_loc_GW+1:div:length(signal);
    signal2=interp1(peak_loc_GW:length(signal),...
        signal(peak_loc_GW:length(signal)),xi2,'linear');

    % Join both interpolations
    signal_interpolated=[signal1 signal2];
    
    % Secure that all signals have same lenght after the interpolation.
    if length(signal_interpolated)<maximum+1
        for j=1:maximum+1-length(signal_interpolated)
            signal_interpolated=[signal_interpolated signal2(length(signal2))];
        end
    end
    
    % Add signal to matrix.
    signals_separates(i,:)=signal_interpolated;
    
end

end

