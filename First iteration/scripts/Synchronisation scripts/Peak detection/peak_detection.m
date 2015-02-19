clear all; close all; clc;

load('acc_shank_RK55.mat');

% Set tuning parameter for peak detection in acceleration signal, that is, 
% threshold and minimum number of samples between standing and walking.
threshold_pos = 1.05;    
threshold_neg = -0.9;
diff_threshold = 0.05; %0.05
gap = 390;

% Find all peaks greater than threshold.
[peak_values_l, peak_locations_l] = findpeaks(a_Z_left_shank_1_C, ...
                                              'minpeakheight', threshold_pos);
[peak_values_r, peak_locations_r] = findpeaks(a_Z_right_shank_1_C, ...
                                              'minpeakheight', threshold_pos);
                                         

plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(1:end-1), diff(a_Z_left_shank_1_C) < diff_threshold, 'r');

figure();
plot(time_GW, a_Z_right_shank_1_C);
hold on;
plot(time_GW(1:end-1), diff(a_Z_right_shank_1_C) < diff_threshold, 'r');

%%
diff_vec_l = diff(a_Z_left_shank_1_C) < diff_threshold;
diff_vec_r = diff(a_Z_right_shank_1_C) < diff_threshold;

ind_vec_l = find(diff_vec_l == 0);
ind_vec_r = find(diff_vec_r == 0);

close all;
stem(diff(ind_vec_l));

figure();
stem(diff(ind_vec_r));

%%
close all;
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(find((a_Z_left_shank_1_C < 1.1)&(a_Z_left_shank_1_C > 0.8))), a_Z_left_shank_1_C(find((a_Z_left_shank_1_C < 1.1)&(a_Z_left_shank_1_C > 0.8))), 'r');
%%

% Detect the last peak before a gap of more than gap samples.
last_first_peaks_l = ind_vec_l(diff(ind_vec_l) > gap);
last_first_peaks_r = ind_vec_r(diff(ind_vec_r) > gap);

close all;

plot(time_GW, a_Z_left_shank_1_C);
hold on;
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(last_first_peaks_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

figure();
plot(time_GW, a_Z_right_shank_1_C);
hold on;
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(last_first_peaks_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

%%

% Store every second peak, that is the one at the end of the first step 
% onto the force plate.
last_peaks_l = last_first_peaks_l(1:2:length(last_first_peaks_l));
last_peaks_r = last_first_peaks_r(1:2:length(last_first_peaks_r));

% Store every second peak, that is the one at the end of the first step 
% onto the force plate.
cycle_end_l = last_first_peaks_l(2:2:length(last_first_peaks_l));
cycle_end_r = last_first_peaks_r(2:2:length(last_first_peaks_r));

% Compute number of cycles.
n_cycles = length(last_peaks_l);

% Initialise variables.
sync_peaks_l = 1:n_cycles;
neg_peaks_l = 1:n_cycles;
start_p_l = 1:n_cycles;
end_p_l = 1:n_cycles;
sync_peaks_r = 1:n_cycles;
neg_peaks_r = 1:n_cycles;
start_p_r = 1:n_cycles;
end_p_r = 1:n_cycles;

close all;
subplot(2, 1, 1);
plot(time_GW, a_Z_left_shank_1_C);
hold on;
plot(time_GW(last_peaks_l), a_Z_left_shank_1_C(last_peaks_l), 'r.', 'markersize', 20);
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(last_peaks_l), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');
     
subplot(2, 1, 2);
plot(time_GW, a_Z_right_shank_1_C);
hold on;
plot(time_GW(last_peaks_r), a_Z_right_shank_1_C(last_peaks_r), 'r.', 'markersize', 20);
% Vertical line at the location of sync_peaks.
hx = graph2d.constantline(time_GW(last_peaks_r), 'LineStyle',':', 'LineWidth', 1 , 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

%%

for cycle = 1:9

%---Left shank---%

% Set interval.
start_p_l(cycle) = last_peaks_l(cycle) - 600;
end_p_l(cycle) = last_peaks_l(cycle) + 100;
                                  
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

% Set interval.
start_p_r(cycle) = last_peaks_r(cycle) - 600;
end_p_r(cycle) = last_peaks_r(cycle) + 100;
                                  
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
plot(time_GW(last_peaks_l), a_Z_left_shank_1_C(last_peaks_l), 'g.', 'markersize', 20);
     
subplot(2, 1, 2);
plot(time_GW, a_Z_right_shank_1_C);
hold on;
plot(time_GW(sync_peaks_r), a_Z_right_shank_1_C(sync_peaks_r), 'r.', 'markersize', 30);
plot(time_GW(neg_peaks_r), a_Z_right_shank_1_C(neg_peaks_r), 'm.', 'markersize', 20);
plot(time_GW(last_peaks_r), a_Z_right_shank_1_C(last_peaks_r), 'g.', 'markersize', 20);
