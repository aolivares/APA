% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------------- APA Analysis --------------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accele,rometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
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
% * Last modification: 06/05/2015

% INFORMATION: This file contains the routine to detect when the second
% step happens, determine the APAs of the FP and GW signals in this case 
% and the correlation between them. 

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;

% Set flags which control the visibility of the figures.
showPlotsCheck = 'yes';
showPlotsAPA = 'yes';
showPlotsCorr = 'yes';

% -------------------------------------------------------------------------
% 1) Select the .mat file and extrat the data form timeseries.
% -------------------------------------------------------------------------

% Selection and load of the synchronised file.
[filename, filepath] = uigetfile('*.mat', ...
    'Select the data file from the patient (.mat)', '../../data/Synchronised/');

load(fullfile(filepath,filename));

% Extract data from timeseries for plot.

% Sum of force of all force plate cells.
force_sum_complete_ts = append(force_sum_ts{1, :});
force_sum_data = force_sum_complete_ts.data;

% Force of right-left and front-back feet (four signals).
force_sensors_complete_ts = append(force_sensors_ts{1, :});
fs_data = force_sensors_complete_ts.data;

% Antereo-Posterior center of pressure.
AP_COP_complete_ts = append(AP_COP_ts{1, :});
AP_COP_data = AP_COP_complete_ts.data;

% Medio-Lateral center of pressure.
ML_COP_complete_ts = append(ML_COP_ts{1, :});
ML_COP_data = ML_COP_complete_ts.data;

% Acceleration of the shanks (X-Z axes and left-right direction).
a_shanks_complete_ts = append(a_shanks{1, :});
a_shanks_data = a_shanks_complete_ts.data;

% Acceleration of the trunk (X, Y, Z axes).
a_trunk_complete_ts = append(a_trunk{1, :});
a_trunk_data = a_trunk_complete_ts.data;

% Angular Velocity of trunk (x, Y, Z axes).
g_trunk_complete_ts = append(g_trunk{1, :});
g_trunk_data = g_trunk_complete_ts.data;

% -------------------------------------------------------------------------
% 2) Determine when the second step occurred.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 2.1) Determine the edges (interval) for the second step.
% -------------------------------------------------------------------------

% Determine when the second step occurred. We extract the part of the signal
% when the patient carried out the step. It happens between the beginning of
% activity period of each cycle and the end of FP data.
% We use the edges of the activity detection in the left shank because
% usually used to be clearer. But the difference between edges of the right
% and left shank is very small.

% We consider the second step starts in the beginning of the second period
% of activity of each cycle.
init_second_step = edges_left(3:4:length(edges_left));

% Determine the time point where there isn't force signal (considerating 
% both feet) for each cycle.
force_sum_data = reshape(force_sum_data(3,1,:), 1, ...
                 length(force_sum_data(3,1,:)));

edges_FP_signal = find(diff(force_sum_data > 100)~=0);
final_second_step_edge = edges_FP_signal(2:2:...
                        length(edges_FP_signal));
final_second_step = force_sum_complete_ts.time(final_second_step_edge);

%--------------------------------------------------------------------------
% 3) Detect APA in trunk signal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 3.1) Determine when the second step (when patient goes down from
% plateform) starts with left or right foot.
%--------------------------------------------------------------------------
% When the patient starts to step with left foot, last value of the cycle
% of ML COP is positive because the step finishes the pressure is located
% in the right foot.
% When the patient starts with right foot, the last value of the ML COP
% cycle will be negative.
AP_COP_data =  reshape(AP_COP_data(3, 1, :), ...
                [1, max(size(AP_COP_data))]);
            
ML_COP_data =  reshape(ML_COP_data(3, 1, :), ...
                [1, max(size(ML_COP_data))]);
            
cycle_start_right = find(ML_COP_data(final_second_step_edge) < 0);


%--------------------------------------------------------------------------
% 3.2) Determine the APA peaks in the Center of Pressure.
%--------------------------------------------------------------------------

% Firstly, we detect the APA peaks in Medio-Lateral direction because we
% only have one peak, so is easier to not confuse with the point of time
% when this happens. This point is exactly the same that in the AP
% direction. Therefore, we can use the result to calculate this after.

% If the patient starts to walk with the right foot, the APA peak is
% positive. However, if the patient starts to walk with the left foot,
% the APA peak is negative.

% Determine the initial and end point (time) of each interval where 
% we need to find the peaks, i.e, the second activity period of each cycle.
initcross = time_GW(init_second_step);
finalcross = final_second_step';

% Calculate the peaks in each interval.
for k = 1:length(initcross)
    
     % To find a noninteger value, we use a tolerance value based on our data.
     % Otherwise, the result is sometimes an empty matrix due to 
     % floating-point roundoff error.
     initcross_COP = find(abs(ML_COP_complete_ts.time - initcross(k))...
                    < 0.001);
     finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k))...
                    < 0.001);
     value_initcross_AP(k) = AP_COP_data(initcross_COP);
     
     % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a positive peak.
        
      % Find the longest positive peak.
      [ peaks_index ] = findMaxPeaks(ML_COP_data, initcross_COP,...
                        finalcross_COP,1);
      peaks_APA_ML_COP(k) =  peaks_index;
      
    else % Look for a negative peak.
        
         % Find longest negative peak. 
          [ peaks_index ] = findMaxPeaks(ML_COP_data, initcross_COP,...
                        finalcross_COP,2);
          peaks_APA_ML_COP(k) =  peaks_index;
                           
    end
    
   % Find longest negative peak around the peak calculated above because 
   % this happens almost at the same time and it is most accurate in the 
   % lateral direction. 
  [ peaks_index ] = findMaxPeaks(AP_COP_data, peaks_APA_ML_COP(k)-30,...
                peaks_APA_ML_COP(k)+30,2);
            
  peaks_APA_AP_COP(k) =  peaks_index;
  
  % Find the next longest positive peak to determine the height of the peak. 
  [ peaks_index ] = findMaxPeaks(AP_COP_data, peaks_APA_AP_COP(k),...
                peaks_APA_AP_COP(k)+100,1);
  peaks_APA_AP_COP_2(k) = peaks_index;
  
  % Find the longest negative peak( this is when the COP is in the stance
  % foot).
  [ peaks_index ] = findMaxPeaks(AP_COP_data, peaks_APA_AP_COP_2(k),...
                finalcross_COP,2);
            
  peaks_APA_AP_COP_3(k) = peaks_index;  
end

% We calculate the value of the APA peaks.
value_APA_ML_COP = ML_COP_data(peaks_APA_ML_COP);
value_APA_AP_COP = AP_COP_data(peaks_APA_AP_COP);
value_APA_AP_COP_2 = AP_COP_data(peaks_APA_AP_COP_2);
value_APA_AP_COP_3 = AP_COP_data(peaks_APA_AP_COP_3);

%------------------------------- Plots-------------------------------------
if strcmpi(showPlotsCheck,'yes')
subplot(2,1,1)
plot(AP_COP_complete_ts.time, AP_COP_data, 'g');
hold on;
plot(AP_COP_complete_ts.time(peaks_APA_AP_COP), value_APA_AP_COP , 'r.');
plot(AP_COP_complete_ts.time(peaks_APA_AP_COP_2), value_APA_AP_COP_2 , 'y.');
plot(AP_COP_complete_ts.time(peaks_APA_AP_COP_3), value_APA_AP_COP_3 , 'b.');
% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['AP COP with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('AP COP (mm)');

subplot(2,1,2)
plot(ML_COP_complete_ts.time, ML_COP_data, 'g');
hold on;
plot(ML_COP_complete_ts.time(peaks_APA_ML_COP), value_APA_ML_COP , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['ML COP with lines marker when the' ...
       'patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('ML COP (mm)');
end

% %--------------------------------------------------------------------------
% % 3.3) Determine the next peak after the APA peaks detected in ML direction. 
% % This point will indicate the interval to find the APA peaks in the 
% % acceleration signals 
% %--------------------------------------------------------------------------
% initcross_COP = peaks_APA_ML_COP;
% 
% % Calculate the peaks in each interval.
% for k = 1:length(initcross)
%     % To find a noninteger value, we use a tolerance value based on our data.
%     % Otherwise, the result is sometimes an empty matrix due to 
%     % floating-point round off error.
%     finalcross_COP = find(abs(ML_COP_complete_ts.time - finalcross(k))...
%                     < 0.001);
%                 
%      % Differenciate when the patient starts with left or right foot.
%     if (find(cycle_start_right == k))% Look for a negative peak.
%         
%          % Find all peaks in each interval.
%         [neg_peak_values, neg_peak_locations] = findpeaks(...
%                                 -ML_COP_data(...
%                                 initcross_COP(k):finalcross_COP));
% 
%         % Store the index of the first negative peak.
%         neg_peak_locations = sort(neg_peak_locations);
%         finalcross_acc(k) = neg_peak_locations(1) + initcross_COP(k) -1;
% 
%                     
%     else % Look for a positive peak.
%         
%          % Find all peaks in each interval.
%         [pos_peak_values, pos_peak_locations] = findpeaks(...
%                                 ML_COP_data(...
%                                 initcross_COP(k):finalcross_COP));
% 
%         % Store the index of the first negative peak.
%         pos_peak_locations = sort(pos_peak_locations);
%         finalcross_acc(k) = pos_peak_locations(1) + initcross_COP(k) -1;
%                            
%     end
% end

                    
%--------------------------------------------------------------------------
% 3.4) Determine the APA peaks in the Acceleration Signals.
%--------------------------------------------------------------------------
% Calculate the indexes to be able to use them for the accelerations
% signals.
% initcross_acc_time = ML_COP_complete_ts.time (initcross_COP);
% finalcross_acc_time = ML_COP_complete_ts.time (finalcross_acc);

a_trunk_data_X = reshape(a_trunk_data(1, 1, :), ...
                [1, max(size(a_trunk_data))]);
a_trunk_data_Y = reshape(a_trunk_data(2, 1, :), ...
                [1, max(size(a_trunk_data))]);
            
% Apply a lowpass filter to the accelerameter signals.
% Definition of the filter's parameters.
Fs = 200;
fc = 2;

% FFT to check.
L =length(a_trunk_data_X);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(a_trunk_data_X,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure()
stem(f,2*abs(Y(1:NFFT/2+1)));

% Desing of FIR filter. The first parameter is the order of the filter. The
% next one represents the cutoff frecuency. This can be a value between 0
% and 1, where 1 is the Nyquist frecuency (sample rate/2).
b=fir1(30,fc/(Fs/2));
figure()
freqz(b);

a_trunk_data_X = filter(b,1,a_trunk_data_X);
a_trunk_data_Y = filter(b,1,a_trunk_data_Y); 

% Show the fft of the signals after the filtering.
Y = fft(a_trunk_data_X,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure()
stem(f,2*abs(Y(1:NFFT/2+1)));


for k = 1:length(initcross)
    
  % We obtain the interval where the APA peaks in the Acc signals appear.
  initcross_acc = find(abs( a_trunk_complete_ts.time- initcross(k))...
                    < 0.001)+20;
  finalcross_acc = find(abs( a_trunk_complete_ts.time- finalcross(k))...
                    < 0.001) - 15;
                
  % We ontain the value of the Acc peak in the ML direction. To do that, we 
  % use a small interval around the limit calculated before (e.i, the second
  % period of the APA).
                
   % Differenciate when the patient starts with left or right foot.
    if (find(cycle_start_right == k))% Look for a positive peak.
        
      % Find the longest positive peak.
      [ peaks_index ] = findMaxPeaks(a_trunk_data_Y, initcross_acc,...
                        finalcross_acc,1);
      peaks_APA_acc_Y(k) =  peaks_index;
      
      % Find the lext negative peak.
      [ peaks_index ] = findMaxPeaks(a_trunk_data_Y, peaks_APA_acc_Y(k),...
                        finalcross_acc,2);
      peaks_APA_acc_Y_2(k) =  peaks_index;     
      
    else % Look for a negative peak.
        
         % Find longest negative peak. 
          [ peaks_index ] = findMaxPeaks(a_trunk_data_Y, initcross_acc,...
                            finalcross_acc,2);
          peaks_APA_acc_Y(k) =  peaks_index;
          
         % Find the next positive peak.
          [ peaks_index ] = findMaxPeaks(a_trunk_data_Y, peaks_APA_acc_Y(k),...
                            finalcross_acc,1);
          peaks_APA_acc_Y_2(k) =  peaks_index;
          
    end
      
  % Find the longest negative peak in AP direction.
%   [ peaks_index ] = findMaxPeaks(a_trunk_data_X, initcross_acc,...
%                         finalcross_acc, 2);
%   peaks_APA_acc_X(k) = peaks_index; 

   [neg_peak_values, neg_peak_locations] = min(...
                                    a_trunk_data_X(initcross_acc:finalcross_acc));
   peaks_APA_acc_X(k) = neg_peak_locations + initcross_acc -1;
  
 end


% We calculate the value of the APA peaks.
value_APA_acc_X = a_trunk_data_X(peaks_APA_acc_X);
value_APA_acc_Y = a_trunk_data_Y(peaks_APA_acc_Y);
value_APA_acc_Y_2 = a_trunk_data_Y(peaks_APA_acc_Y_2);
%------------------------------- Plots figure------------------------------
if strcmpi(showPlotsCheck,'yes')
figure ()
subplot(2,1,1)
plot(a_trunk_complete_ts.time, a_trunk_data_X, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_acc_X), value_APA_acc_X , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the X-axis of the trunk with lines marker when the' ...
       ' patient steps with the second time and the APA peaks']);
xlabel('Time in s');
ylabel('Acceleration (g)');

subplot(2,1,2)
plot(a_trunk_complete_ts.time, a_trunk_data_Y, 'g');
hold on;
plot(a_trunk_complete_ts.time(peaks_APA_acc_Y), value_APA_acc_Y , 'r.');
plot(a_trunk_complete_ts.time(peaks_APA_acc_Y_2), value_APA_acc_Y_2 , 'y.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Acceleration of the Y-axis of the trunk with lines marker when the' ...
       ' patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Acceleration in (g)');
end

%--------------------------------------------------------------------------
% 3.5) Determine the APA peaks in the Gyroscope Signals.
%--------------------------------------------------------------------------
g_trunk_data_X = reshape(g_trunk_data(1, 1, :), ...
                [1, max(size(g_trunk_data))]);
g_trunk_data_Y = reshape(g_trunk_data(2, 1, :), ...
                [1, max(size(g_trunk_data))]);
            
%------------------------------- Plots-------------------------------------
% Plot to check the gyroscope signal
if strcmpi(showPlotsCheck,'yes')
figure ()
subplot(2,1,1)
plot(g_trunk_complete_ts.time, g_trunk_data_X, 'g');
% hold on;
% plot(a_trunk_complete_ts.time(peaks_APA_acc_X), value_APA_acc_X , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Angular Velocity of the X-axis of the trunk with lines marker when the' ...
       ' patient steps with the second time and the APA peaks']);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');
axis([170, 210, -100, 100]);

subplot(2,1,2)
plot(g_trunk_complete_ts.time, g_trunk_data_Y, 'g');
% hold on;
% plot(a_trunk_complete_ts.time(peaks_APA_acc_Y), value_APA_acc_Y , 'r.');

% Vertical line init.
hx = graph2d.constantline(time_GW(init_second_step), 'LineStyle',':',...
    'LineWidth', 2 , 'Color', 'r');
changedependvar(hx,'x');

% Vertical line final.
hx = graph2d.constantline(final_second_step, 'LineStyle',':', ...
    'LineWidth', 2 , 'Color', 'm');
changedependvar(hx,'x');

title(['Angular Velocity of the Y-axis of the trunk with lines marker when the' ...
       ' patient steps with the second time and the APA peak' ]);
xlabel('Time in s');
ylabel('Angular Velocity (º/s)');
axis([170, 210, -40, 40]);
end

 %-------------------------------------------------------------------------
 % 4) Correlations
 %-------------------------------------------------------------------------
 
 % Correlation in AP direction.
 value_APA_AP_COP_1 = abs(value_APA_AP_COP - value_initcross_AP);
 value_APA_acc_X_1 = abs(value_APA_acc_X- mode(a_trunk_data_X));
 
 [corr_trunk_AP_1, prob_trunk_AP_1] = corr(value_APA_acc_X_1',...
                                    value_APA_AP_COP_1');

 value_APA_AP_COP_2 = abs( value_APA_AP_COP - value_APA_AP_COP_2);
 [corr_trunk_AP_2, prob_trunk_AP_2] = corr(value_APA_acc_X_1',...
                                    value_APA_AP_COP_2');
% Correlation in ML direction.
value_APA_ML_COP_1 = abs(value_APA_ML_COP);
value_APA_acc_Y_1 = abs(value_APA_acc_Y) - mode(a_trunk_data_Y);

[corr_trunk_ML_1, prob_trunk_ML_1] = corr(value_APA_acc_Y_1',...
                                    value_APA_ML_COP_1');
 
% -------------------------- Correlations Plots----------------------------
if strcmpi(showPlotsCorr,'yes')
figure ()

plot(value_APA_acc_X_1, value_APA_AP_COP_2, '.r');
title('Linear correlation between peak AP_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('AP_COP (mmm)');

figure ()

plot(value_APA_acc_Y_1, value_APA_ML_COP_1, '.r');
title('Linear correlation between peak ML_COP and peak Acc trunk');
xlabel('Acceleration (g)');
ylabel('AP_COP (mmm)');

end

%-------------------------------------------------------------------------
% 5) Trajectory during a cycle of APA.
%-------------------------------------------------------------------------                               
 
 %-------------------------------------------------------------------------
 % 5.1) Center of Pressure (COP) Trajectory.
 %-------------------------------------------------------------------------
 % We are goint to show the trajectory of the APA in the seventh cycle.
   initcross_7 = find(abs(AP_COP_complete_ts.time - initcross(7))...
                    < 0.001);
   finalcross_7 = find(abs(AP_COP_complete_ts.time - finalcross(7))...
                    < 0.001);
                
 %--------------------------- Plot figure----------------------------------
  if strcmpi(showPlotsAPA,'yes')
  figure()
  plot3(ML_COP_data(initcross_7:finalcross_7),...
                     AP_COP_data(initcross_7:finalcross_7),...
                     AP_COP_complete_ts.time(initcross_7:finalcross_7));
                 
  title('Trajectory of the COP during APA (7th cycle) ' );
  xlabel('ML COP (mm)');
  ylabel('AP COP (mm)');
  zlabel('Time (s)');
  grid on
  axis square
  end
   % We are going to show the trajectory of the APA in the eithth cycle.
   initcross_8 = find(abs(AP_COP_complete_ts.time - initcross(8))...
                    < 0.001);
   finalcross_8 = find(abs(AP_COP_complete_ts.time - finalcross(8))...
                    < 0.001);
                
 % --------------------------- Plot figure --------------------------------
  if strcmpi(showPlotsAPA,'yes')
  figure()
  plot3(ML_COP_data(initcross_8:finalcross_8),...
                     AP_COP_data(initcross_8:finalcross_8),...
                     AP_COP_complete_ts.time(initcross_8:finalcross_8));
                 
  title('Trajectory of the COP during APA (8th cycle) ' );
  xlabel('ML COP (mm)');
  ylabel('AP COP (mm)');
  zlabel('Time (s)');
  grid on
  axis square
  end
  
 %-------------------------------------------------------------------------
 % 5.2) Acceleration (Acc) Trajectory.
 %-------------------------------------------------------------------------
 % We are goint to show the trajectory of the APA in the second cycle.
  initcross_7 = find(abs(a_trunk_complete_ts.time - initcross(7))...
                    < 0.001);
  finalcross_7 = find(abs(a_trunk_complete_ts.time - finalcross(7))...
                    < 0.001);
                
  % -------------------------- Plot figure---------------------------------
  if strcmpi(showPlotsAPA,'yes')
  figure()
  plot3(a_trunk_data_Y(initcross_7:finalcross_7),...
                     a_trunk_data_X(initcross_7:finalcross_7),...
                     a_trunk_complete_ts.time(initcross_7:finalcross_7));
                 
  title('Trajectory of the Acceleration during APA ' );
  xlabel('Acc Y (g)');
  ylabel('Acc X (g)');
  zlabel('Time (s)');
  grid on
  axis square
  end
  