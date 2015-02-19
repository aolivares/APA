clear all; close all; clc;

load('acc_shank_RK55.mat');

% Set tuning parameter for peak detection in acceleration signal, that is, 
% threshold and minimum number of samples between standing and walking.
threshold_pos = 1.05;
threshold_neg = -0.9;
gap = 500;

% Find all peaks greater than threshold.
[peak_values_l, peak_locations_l] = findpeaks(a_Z_left_shank_1_C, ...
                                              'minpeakheight', threshold_pos);
[peak_values_r, peak_locations_r] = findpeaks(a_Z_right_shank_1_C, ...
                                              'minpeakheight', threshold_pos);

% Detect the last peak before a gap of more than gap samples.
last_first_peaks_l = peak_locations_l(diff(peak_locations_l) > gap);
last_first_peaks_r = peak_locations_r(diff(peak_locations_r) > gap);

% Store every second peak, that is the one at the end of the first step 
% onto the force plate.
last_peaks_l = last_first_peaks_l(1:2:length(last_first_peaks_l));
last_peaks_r = last_first_peaks_r(1:2:length(last_first_peaks_r));

% Compute number of cycles.
n_cycles = length(last_peaks_l);

% Initialise variables.
sync_peaks_l = 1:n_cycles;
neg_peaks_l = 1:n_cycles;
start_p_l = 1:n_cycles;
end_p_l = 1:n_cycles;

for cycle = 1:n_cycles

start_p_l(cycle) = last_peaks_l(cycle) - 800;
end_p_l(cycle) = last_peaks_l(cycle) + 400;
                                  
% Find all peaks smaller than threshold in the interval specified by start_p and end_p.
[neg_peak_values_int_l, neg_peak_locations_int_l] = findpeaks(-a_Z_left_shank_1_C(start_p_l(cycle):end_p_l(cycle)), ...
                                              'minpeakheight', threshold_neg);

% Store the index of the highest peak in a.                                      
neg_peaks_l(cycle) = find(a_Z_left_shank_1_C(start_p_l(cycle):end_p_l(cycle)) == -max(neg_peak_values_int_l)) + start_p_l(cycle) - 1;
                                     
 
 % Find all peaks greater than threshold in the interval specified by a and
 % end_p
[peak_values_int_l, peak_locations_int_l] = findpeaks(a_Z_left_shank_1_C(neg_peaks_l(cycle):end_p_l(cycle)), ...
                                              'minpeakheight', threshold_pos);

% Store the index of the highest peak in a.                                      
sync_peaks_l(cycle) = find(a_Z_left_shank_1_C(neg_peaks_l(cycle):end_p_l(cycle)) == max(peak_values_int_l), 1, 'first') + neg_peaks_l(cycle) - 1;

end


close all;

plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(sync_peaks_l), a_Z_left_shank_1_C(sync_peaks_l), 'r.', 'markersize', 40);
plot(time_GW(neg_peaks_l), a_Z_left_shank_1_C(neg_peaks_l), 'm.', 'markersize', 20);
plot(time_GW(last_peaks_l), a_Z_left_shank_1_C(last_peaks_l), 'g.', 'markersize', 20);
     

