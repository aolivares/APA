clear all; close all; clc;

load('acc_shank_ES39.mat');

% Set tuning parameter for peak detection in acceleration signal, that is, 
% threshold and minimum number of samples between standing and walking.
threshold_pos = 1.05;    
threshold_neg = -0.9;
gap = 1313;

% Find all peaks greater than threshold.
[peak_values_l, peak_locations_l] = findpeaks(a_Z_left_shank_1_C, ...
                                              'minpeakheight', threshold_pos);
[peak_values_r, peak_locations_r] = findpeaks(a_Z_right_shank_1_C, ...
                                              'minpeakheight', threshold_pos);
          
% Compute distance between two peaks.                             
peak_distance_l = diff(peak_locations_l);
peak_distance_r = diff(peak_locations_r);


%%


close all;
stem(peak_distance_l);

figure();
stem(peak_distance_r);

%%

% Create logical vector to select last peak of each cycle. That is the one
% before a gap of at least the number of samples stored in gap,
% respectively. 
select_last_l = peak_distance_l > gap;
select_last_r = peak_distance_r > gap;

last_peaks_l = peak_locations_l(select_last_l);
last_peaks_r = peak_locations_r(select_last_r);

% Shift logical vector by one to the right to select the first instead of
% the last peak of a cycle.
select_first_l = [0, select_last_l(1:length(select_last_l) - 1)];
select_first_r = [0, select_last_r(1:length(select_last_r) - 1)];


start_p_l  = [peak_locations_l(1), peak_locations_l(logical(select_first_l))]-100;
start_p_r  = [peak_locations_r(1), peak_locations_r(logical(select_first_r))]-100;

end_p_l  = [peak_locations_l(1), peak_locations_l(logical(select_first_l))]+400;
end_p_r  = [peak_locations_r(1), peak_locations_r(logical(select_first_r))]+400;
                                          

plot(time_GW, a_Z_left_shank_1_C);
hold on;
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(start_p_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(end_p_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

figure();
plot(time_GW, a_Z_right_shank_1_C);
hold on;
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(start_p_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(end_p_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');


% Compute number of cycles.
n_cycles = length(last_peaks_l);


% Initialise variables.
sync_peaks_l = 1:n_cycles;
neg_peaks_l = 1:n_cycles;
sync_peaks_r = 1:n_cycles;
neg_peaks_r = 1:n_cycles;




%%

for cycle = 1:n_cycles

%---Left shank---%

% Find all peaks smaller than threshold in the interval specified by start_p and end_p.
[neg_peak_values_int_l, neg_peak_locations_int_l] = findpeaks(-a_Z_left_shank_1_C(start_p_l(cycle):end_p_l(cycle)), ...
                                              'minpeakheight', threshold_neg);

% Store the index of the highest peak in a.                                      
neg_peaks_l(cycle) = find(a_Z_left_shank_1_C(start_p_l(cycle):end_p_l(cycle)) == -max(neg_peak_values_int_l), 1, 'last') + start_p_l(cycle) - 1;
                                     
 
 % Find all peaks greater than threshold in the interval specified by a and
 % end_p
[peak_values_int_l, peak_locations_int_l] = findpeaks(a_Z_left_shank_1_C(neg_peaks_l(cycle):end_p_l(cycle)), ...
                                              'minpeakheight', threshold_pos);

% Store the index of the highest peak in a.                                      
sync_peaks_l(cycle) = find(a_Z_left_shank_1_C(neg_peaks_l(cycle):end_p_l(cycle)) == max(peak_values_int_l), 1, 'first') + neg_peaks_l(cycle) - 1;


%---Right shank---%
                                  
% Find all peaks smaller than threshold in the interval specified by start_p and end_p.
[neg_peak_values_int_r, neg_peak_locations_int_r] = findpeaks(-a_Z_right_shank_1_C(start_p_r(cycle):end_p_r(cycle)), ...
                                              'minpeakheight', threshold_neg);

% Store the index of the highest peak in a.                                      
neg_peaks_r(cycle) = find(a_Z_right_shank_1_C(start_p_r(cycle):end_p_r(cycle)) == -max(neg_peak_values_int_r), 1, 'last') + start_p_r(cycle) - 1;
                                     
 
 % Find all peaks greater than threshold in the interval specified by a and
 % end_p
[peak_values_int_r, peak_locations_int_r] = findpeaks(a_Z_right_shank_1_C(neg_peaks_r(cycle):end_p_r(cycle)), ...
                                              'minpeakheight', threshold_pos);

% Store the index of the highest peak in a.                                      
sync_peaks_r(cycle) = find(a_Z_right_shank_1_C(neg_peaks_r(cycle):end_p_r(cycle)) == max(peak_values_int_r), 1, 'first') + neg_peaks_r(cycle) - 1;


end


close all;

subplot(2, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l), a_Z_left_shank_1_C(sync_peaks_l), 'r.', 'markersize', 30);
plot(time_GW(neg_peaks_l), a_Z_left_shank_1_C(neg_peaks_l), 'm.', 'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(start_p_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(end_p_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
     
subplot(2, 1, 2);
plot(time_GW, a_Z_right_shank_1_C);
hold on;
plot(time_GW(sync_peaks_r), a_Z_right_shank_1_C(sync_peaks_r), 'r.', 'markersize', 30);
plot(time_GW(neg_peaks_r), a_Z_right_shank_1_C(neg_peaks_r), 'm.', 'markersize', 20);

% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(start_p_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(end_p_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
