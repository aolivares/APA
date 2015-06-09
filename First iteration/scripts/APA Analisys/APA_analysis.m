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
% * Last modification: 08/06/2015

% INFORMATION: This file contains the routine to detect when the second
% step happens, determine the APAs of the FP and GW signals in this case 
% and the correlation between them. 

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;
gw = gwLibrary;

% Set extra-caculations
peakManually = 'no';
PCA = 'yes';

% -------------------------------------------------------------------------
% 1) Select, read and obtain information from the excel file.
% -------------------------------------------------------------------------
 
% Select only one data files with a dialog box (only .xlsl files).
[filename_excel, filepath] = uigetfile('*.xlsx', ...
    'Select the Excel file with the data (.xlsx)');

% Read data from *.xlsx where are stored all filenames and other 
% interesting information.
[~,file_excel] = xlsread([filepath, filename_excel]);
[rows,columns] = size(file_excel);

% The names of the files are in the first column and draw from the forth row.
for i = 4:rows
    
    % Convert the namefile in a char type to be able to read the file.
    filename = char(file_excel(i));
    
    % Data load.
    load(fullfile('../../data/Synchronised/',filename));
    
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
    
    % Determine the initial and end point (time) of each interval where 
    % we need to find the peaks, i.e, the second activity period of each cycle.
    initcross = time_GW(init_second_step);
    finalcross = final_second_step';
    
    % Determine the initial and end point of each interval for FP and GW.
    for k = 1:length(initcross)
    
     % To find a noninteger value, we use a tolerance value based on our data.
     % Otherwise, the result is sometimes an empty matrix due to 
     % floating-point roundoff error.
     initcross_COP (k) = find(abs(ML_COP_complete_ts.time - initcross(k))...
                    < 0.001);
     finalcross_COP(k) = find(abs(ML_COP_complete_ts.time - finalcross(k))...
                    < 0.001);
                
     initcross_acc (k)= find(abs( a_trunk_complete_ts.time- initcross(k))...
                    < 0.001);
     finalcross_acc (k)= find(abs( a_trunk_complete_ts.time- finalcross(k))...
                    < 0.001);
    end
    
    % Obtaint he signals from GW.
    a_trunk_data_X = reshape(a_trunk_data(1, 1, :), ...
                [1, max(size(a_trunk_data))]);
    a_trunk_data_Y = reshape(a_trunk_data(2, 1, :), ...
                    [1, max(size(a_trunk_data))]);

    g_trunk_data_X = reshape(g_trunk_data(1, 1, :), ...
                    [1, max(size(g_trunk_data))]);
    g_trunk_data_Y = reshape(g_trunk_data(2, 1, :), ...
                    [1, max(size(g_trunk_data))]);
                
%--------------------------------------------------------------------------
% 3.1) Lowpass Filter.
%--------------------------------------------------------------------------
    % Apply a lowpass filter to the accelerameter  and gyroscope signals.
    % Definition of the filter's parameters.
    Fs = 200;
    fc = 2;

    % FFT to check.
    L =length(a_trunk_data_X);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(a_trunk_data_X, NFFT)/L;
    f = Fs/2*linspace(0, 1, NFFT/2 + 1);

    % Desing of FIR filter. The first parameter is the order of the filter. The
    % next one represents the cutoff frecuency. This can be a value between 0
    % and 1, where 1 is the Nyquist frecuency (sample rate/2).
    b=fir1(30, fc/(Fs/2));

    a_trunk_data_X = filter(b,1,a_trunk_data_X);
    a_trunk_data_Y = filter(b,1,a_trunk_data_Y); 

    g_trunk_data_X = filter(b,1,g_trunk_data_X);
    g_trunk_data_Y = filter(b,1,g_trunk_data_Y);
    
% -------------------------------------------------------------------------
% 3.3) Synchronisation of all cycles.
% -------------------------------------------------------------------------

    % The goal of this part of the code is to carry out the mean of all cycles
    % where the patient repeated the same protocol. To do this, we use the cross
    % correlation when the patient does the second step.

    % First interation to align the signals and carry out the mean between both.
    % If patient starts with left foot, we invert the signal to carry out
    % the mean.
    if (find(cycle_start_right == 1))sing1=1; else sing1=-1; end
    if (find(cycle_start_right == 2))sing2=1; else sing2=-1; end

    % Center of pressure   
    ML_COP_mean = aligned_signals( sing1.*ML_COP_data(initcross_COP(2):finalcross_COP(2)),...
                         sing2.*ML_COP_data(initcross_COP(1):finalcross_COP(1)));
    AP_COP_mean = aligned_signals( AP_COP_data(initcross_COP(2):finalcross_COP(2)),...
                         AP_COP_data(initcross_COP(1):finalcross_COP(1)));

    % Acceleration.
    a_trunk_data_X_mean = aligned_signals( a_trunk_data_X(initcross_acc(1):finalcross_acc(1)),...
                         a_trunk_data_X(initcross_acc(2):finalcross_acc(2)));

    a_trunk_data_Y_mean = aligned_signals(sing1.* a_trunk_data_Y(initcross_acc(1):finalcross_acc(1)),...
                         sing2.*a_trunk_data_Y(initcross_acc(2):finalcross_acc(2)));

    % Angular Velocity.
    g_trunk_data_X_mean = aligned_signals( g_trunk_data_X(initcross_acc(1):finalcross_acc(1)),...
                         g_trunk_data_X(initcross_acc(2):finalcross_acc(2)));

    g_trunk_data_Y_mean = aligned_signals( sing1.*g_trunk_data_Y(initcross_acc(1):finalcross_acc(1)),...
                         sing2.*g_trunk_data_Y(initcross_acc(2):finalcross_acc(2)));


     for k = 3:length(initcross) 

        % If patient starts with left foot, we invert the signal to carry out
        % the mean.
        if (find(cycle_start_right == k))sing=1;else sing=-1;end

        % Align the signals and carry out the mean.
        % We consider the sing in ML direction.
         ML_COP_mean = aligned_signals(sing* ML_COP_data(initcross_COP(k):finalcross_COP(k)),...
                    ML_COP_mean );

         a_trunk_data_Y_mean = aligned_signals(sing* a_trunk_data_Y(initcross_acc(k):finalcross_acc(k)),...
                            a_trunk_data_Y_mean );

         g_trunk_data_Y_mean = aligned_signals( sing*g_trunk_data_Y(initcross_acc(k):finalcross_acc(k)),...
                            g_trunk_data_Y_mean );


         % In AP direcction, there isn't change of sing.
          AP_COP_mean = aligned_signals( AP_COP_data(initcross_COP(k):finalcross_COP(k)),...
                             AP_COP_mean);

          a_trunk_data_X_mean = aligned_signals( a_trunk_data_X(initcross_acc(k):finalcross_acc(k)),...
                             a_trunk_data_X_mean);

          g_trunk_data_X_mean = aligned_signals( g_trunk_data_X(initcross_acc(k):finalcross_acc(k)),...
                             g_trunk_data_X_mean);

     end


    % Select the APA points in all signals. You can do this manually or
    % automatically.
     if strcmpi(peakManually,'yes')
    index_APA_AP_COP = gw.getDCindexes(AP_COP_mean,'Select the three APA points (negative-positive-negative)');
    close(gcf)

    index_APA_ML_COP = gw.getDCindexes(ML_COP_mean,'Select the APA point(logest positive)');
    close(gcf)

    index_APA_acc_X = gw.getDCindexes(a_trunk_data_X_mean,'Select the APA point(logest negative)');
    close(gcf)

    index_APA_acc_Y = gw.getDCindexes(a_trunk_data_Y_mean,'Select the two APA points(positive-negative)');
    close(gcf)

    index_APA_gyro_X = gw.getDCindexes(g_trunk_data_X_mean,'Select the APA point(longest negative)');
    close(gcf)

    index_APA_gyro_Y = gw.getDCindexes(g_trunk_data_Y_mean,'Select the two APA points (positive-negative)');
    close(gcf)

    % Calculate the APA Parameters.
    APA_COP_AP_1 = abs( AP_COP_mean(index_APA_AP_COP(1) ) - AP_COP_mean(1));
    APA_COP_AP_2 = abs(AP_COP_mean(index_APA_AP_COP(1) ) - AP_COP_mean(index_APA_AP_COP(2) ));
    APA_COP_AP_3 = abs(AP_COP_mean(index_APA_AP_COP(2) ) - AP_COP_mean(index_APA_AP_COP(3) ));

    APA_COP_ML_1 = ML_COP_mean(index_APA_ML_COP);

    APA_Acc_X_1 =  abs( a_trunk_data_X_mean(index_APA_acc_X) - a_trunk_data_X_mean(1));

    APA_Acc_Y_1 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y(1)) - a_trunk_data_Y_mean(1));
    APA_Acc_Y_2 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y(1)) - a_trunk_data_Y_mean(index_APA_acc_Y(2)));
    APA_Acc_Y_3 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y(2)) - a_trunk_data_Y_mean(1));

    APA_Gyro_X_1 =  abs( g_trunk_data_X_mean(index_APA_gyro_X) - g_trunk_data_X_mean(1));

    APA_Gyro_Y_1 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y(1)) - g_trunk_data_Y_mean(1));
    APA_Gyro_Y_2 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y(2)) - g_trunk_data_Y_mean(index_APA_gyro_Y(1)));
    APA_Gyro_Y_3 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y(2)) - g_trunk_data_Y_mean(1));

    APA_COP_duration = abs(AP_COP_complete_ts.time(index_APA_ML_COP) - ...
        AP_COP_complete_ts.time(initcross_COP(1)+length(ML_COP_mean)-1));

    APA_Acc_duration = a_trunk_complete_ts.time(index_APA_acc_Y(2)) - ...
        a_trunk_complete_ts.time(index_APA_acc_Y(1));

    APA_Gyro_duration = a_trunk_complete_ts.time(index_APA_gyro_Y(2)) - ...
        a_trunk_complete_ts.time(index_APA_gyro_Y(1));

    APA_Parameters = [APA_COP_AP_1, APA_COP_AP_2, APA_COP_AP_3, APA_COP_ML_1,...
                     APA_Acc_X_1, APA_Acc_Y_1, APA_Acc_Y_2, APA_Acc_Y_3,...
                     APA_Gyro_X_1, APA_Gyro_Y_1, APA_Gyro_Y_2, APA_Gyro_Y_3,...
                     APA_COP_duration, APA_Acc_duration, APA_Gyro_duration];
                 
     else % The APA peaks are detected automatically.
     
    % ML COP signal.
    % Detect when there is a strong change of level in the signal.
     d=diff(ML_COP_mean);
     [~, neg_peak_locations] = min(d(1:length(d)-50));
     
    % Obtain the APA peak in ML direction. This is the maximum value before
    % this, close to the descent.
     [~, pos_peak_locations] = max(ML_COP_mean(neg_peak_locations-50:neg_peak_locations));
     index_APA_ML_COP_1 = pos_peak_locations + neg_peak_locations - 50;
     
    % AP COP signal.
    % Detect when there is a strong change of level in the signal. We have
    % to detect two changes here.
    d=diff(AP_COP_mean);
    [~, neg_peak_locations] = min(d(1:length(d)-50));
    [~, pos_peak_locations] = max(d(1:neg_peak_locations));
    
    % One of the APA peak is between the both crossover calculated above.
    % The second one is the minimum value after the descent.
    [~,index_APA_AP_COP_1 ]= max(AP_COP_mean(pos_peak_locations: neg_peak_locations));
    index_APA_AP_COP_1 = index_APA_AP_COP_1+ pos_peak_locations;
    
    [~,index_APA_AP_COP_2] = min(AP_COP_mean(neg_peak_locations:length(AP_COP_mean)));
    index_APA_AP_COP_2 = index_APA_AP_COP_2 + neg_peak_locations;
    
    % ACC and GYRO X signal.
    % We have to find the minumum of this signals. This is when the patient
    % goes forward.
    [~,index_APA_Acc_X]= min(a_trunk_data_X_mean);
    [~,index_APA_Gyro_X]= min(g_trunk_data_X_mean);
    
    % ACC and GYRO Y signal.
    % We have to find the minumum of this signals and after that the prior
    % maximum. This is when the patient moves toward left and right.
    [~,index_APA_Acc_Y_1]= min(a_trunk_data_Y_mean);
    [~,index_APA_Gyro_Y_1]= min(g_trunk_data_Y_mean); 
    
    [~,index_APA_Acc_Y_2 ]= max(a_trunk_data_Y_mean(1:index_APA_Acc_Y_1 ));
    [~,index_APA_Gyro_Y_2 ]= max(g_trunk_data_Y_mean(1:index_APA_Gyro_Y_1 ));
     end

%     figure()
%     plot(g_trunk_complete_ts.time(1:length(g_trunk_data_X_mean)), g_trunk_data_X_mean, 'g');
%     hold on;
%     plot(a_trunk_complete_ts.time(peaks_APA_gyro_X), value_APA_gyro_X , 'r.');
%     plot(AP_COP_complete_ts.time(1:length(AP_COP_mean)), AP_COP_mean, 'g');
%     hold on;
%     plot(AP_COP_complete_ts.time(index_APA_AP_COP_1), AP_COP_mean(index_APA_AP_COP_1) , 'r.');
%     plot(AP_COP_complete_ts.time(index_APA_AP_COP_2), AP_COP_mean(index_APA_AP_COP_2) , 'y.');
end
