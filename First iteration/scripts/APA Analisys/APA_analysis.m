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
% step happens, determine the APAs of the FP and GW signals and the 
% correlation between them. 
%
% * 1) Obtain the data.
% 
% * 2) Determine when the second step occurred.
% 
% * 3) Detect APA in trunk signal.
%
% * 4) Apply PCA.
%
% * 5) Correlations.
%
% * 6) PCA between patients
%
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;
gw = gwLibrary;

% Set extra-caculations.
showPlotsCorr = 'yes';
peakManually = 'no';
PCA = 'yes';

% -------------------------------------------------------------------------
% 1) Select, read and obtain information from the excel file.
% -------------------------------------------------------------------------
 
% Select only one data files with a dialog box (only .xlsl files).
% 'SynchronisedData'.
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
% -------------------------------------------------------------------------
% 3.4) Extract features to characterise APAs.
% -------------------------------------------------------------------------

    % Select the APA points in all signals. You can do this manually or
    % automatically.
    if strcmpi(peakManually,'yes')
    index_APA_AP_COP = gw.getDCindexes(AP_COP_mean,'Select the thwo APA points in AP-COP(positive-negative)');
    close(gcf)

    index_APA_ML_COP = gw.getDCindexes(ML_COP_mean,'Select the APA point in ML-COP(logest positive)');
    close(gcf)

    index_APA_acc_X = gw.getDCindexes(a_trunk_data_X_mean,'Select the APA point in Acc-X(logest negative)');
    close(gcf)

    index_APA_acc_Y = gw.getDCindexes(a_trunk_data_Y_mean,'Select the two APA points in Acc-Y(positive-negative)');
    close(gcf)

    index_APA_gyro_X = gw.getDCindexes(g_trunk_data_X_mean,'Select the APA point in Gyro-X(longest negative)');
    close(gcf)

    index_APA_gyro_Y = gw.getDCindexes(g_trunk_data_Y_mean,'Select the two APA points in Gyro-Y(positive-negative)');
    close(gcf)

    % Calculate the APA Parameters.
    APA_COP_AP_1 = abs( AP_COP_mean(index_APA_AP_COP(1) ) - AP_COP_mean(1));
    APA_COP_AP_2 = abs(AP_COP_mean(index_APA_AP_COP(1) ) - AP_COP_mean(index_APA_AP_COP(2) ));
    APA_COP_AP_3 = abs(AP_COP_mean(index_APA_AP_COP(2) ) - AP_COP_mean(1) );

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

    APA_Parameters (i-3,:) = [APA_COP_AP_1, APA_COP_AP_2, APA_COP_AP_3, APA_COP_ML_1,...
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
    [~,index_APA_acc_X]= min(a_trunk_data_X_mean);
    [~,index_APA_gyro_X]= min(g_trunk_data_X_mean);
    
    % ACC and GYRO Y signal.
    % We have to find the minumum of this signals and after that the prior
    % maximum. This is when the patient moves toward left and right.
    [~,index_APA_acc_Y_1]= min(a_trunk_data_Y_mean);
    [~,index_APA_gyro_Y_1]= min(g_trunk_data_Y_mean); 
    
    [~,index_APA_acc_Y_2 ]= max(a_trunk_data_Y_mean(1:index_APA_acc_Y_1 ));
    [~,index_APA_gyro_Y_2 ]= max(g_trunk_data_Y_mean(1:index_APA_gyro_Y_1 ));
    
    % Calculate the APA Parameters.
    APA_COP_AP_1 = abs( AP_COP_mean(index_APA_AP_COP_1) - AP_COP_mean(1));
    APA_COP_AP_2 = abs(AP_COP_mean(index_APA_AP_COP_2) - AP_COP_mean(index_APA_AP_COP_1));
    APA_COP_AP_3 = abs(AP_COP_mean(index_APA_AP_COP_2) - AP_COP_mean(1));

    APA_COP_ML_1 = ML_COP_mean(index_APA_ML_COP_1);

    APA_Acc_X_1 =  abs( a_trunk_data_X_mean(index_APA_acc_X) - a_trunk_data_X_mean(1));

    APA_Acc_Y_1 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y_1) - a_trunk_data_Y_mean(1));
    APA_Acc_Y_2 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y_1) - a_trunk_data_Y_mean(index_APA_acc_Y_2));
    APA_Acc_Y_3 =  abs( a_trunk_data_Y_mean(index_APA_acc_Y_2) - a_trunk_data_Y_mean(1));

    APA_Gyro_X_1 =  abs( g_trunk_data_X_mean(index_APA_gyro_X) - g_trunk_data_X_mean(1));

    APA_Gyro_Y_1 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y_1) - g_trunk_data_Y_mean(1));
    APA_Gyro_Y_2 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y_2) - g_trunk_data_Y_mean(index_APA_gyro_Y_1));
    APA_Gyro_Y_3 =  abs( g_trunk_data_Y_mean(index_APA_gyro_Y_2) - g_trunk_data_Y_mean(1));

    APA_COP_duration = abs(AP_COP_complete_ts.time(index_APA_ML_COP_1) - ...
        AP_COP_complete_ts.time(initcross_COP(1)+length(ML_COP_mean)-1));

    APA_Acc_duration = a_trunk_complete_ts.time(index_APA_acc_Y_2) - ...
        a_trunk_complete_ts.time(index_APA_acc_Y_1);

    APA_Gyro_duration = a_trunk_complete_ts.time(index_APA_gyro_Y_2) - ...
        a_trunk_complete_ts.time(index_APA_gyro_Y_1);

    APA_Parameters(i-3,:) = [APA_COP_AP_1, APA_COP_AP_2, APA_COP_AP_3, APA_COP_ML_1,...
                     APA_Acc_X_1, APA_Acc_Y_1, APA_Acc_Y_2, APA_Acc_Y_3,...
                     APA_Gyro_X_1, APA_Gyro_Y_1, APA_Gyro_Y_2, APA_Gyro_Y_3,...
                     APA_COP_duration, APA_Acc_duration, APA_Gyro_duration];
    end
%--------------------------------------------------------------------------
% 4) PCA.
%--------------------------------------------------------------------------
if strcmpi(PCA,'yes')
    
% -------------------------------------------------------------------------
% 4.1) PCA in Antero-Posterior Direction.
% -------------------------------------------------------------------------
    % Interpolation of the signals.
    max_length = max([length(AP_COP_mean),length(a_trunk_data_X_mean),...
                length(g_trunk_data_X_mean)]);

    AP_COP_mean = interp1([1:length(AP_COP_mean)],AP_COP_mean,[1:max_length]);
    a_trunk_data_X_mean = interp1([1:length(a_trunk_data_X_mean)],...
                            a_trunk_data_X_mean,[1:max_length]);
    g_trunk_data_X_mean = interp1([1:length(g_trunk_data_X_mean)],...
                            g_trunk_data_X_mean,[1:max_length]);
    % Center the data.
    AP_COP_mean_c = AP_COP_mean - mean(find(AP_COP_mean~=NaN));
    a_trunk_data_X_mean_c = a_trunk_data_X_mean - mean(find(a_trunk_data_X_mean~=NaN));
    g_trunk_data_X_mean_c = g_trunk_data_X_mean - mean(find(g_trunk_data_X_mean~=NaN));

    % Apply PCA.
    X_x = [AP_COP_mean_c; a_trunk_data_X_mean_c; g_trunk_data_X_mean_c]';
    
    % Using 'primcomp' funtion, we can obtain the projection of PCA. 
    % COEFF are the eighvectors, SCORE the projections in the ortogonal
    % space, latent are the eighvalues and tsquare is a measure of
    % probability.
    [COEFF_x,SCORE_x,latent,tsquare] = princomp(X_x,'econ');
    
    % The first componets of the SCORE are the most variability, so we will
    % use them to extract features with more accuracity.
    comp_X_1 = SCORE_x(:,1);
    comp_X_2 = SCORE_x(:,2);
    
    % Feature extration of the componentes in AP direction.
    % First component.
    % Detect when there is a strong change of level in the signal. We have
    % to detect two changes here.
    d=diff(comp_X_1 );
    [~, neg_peak_locations] = min(d(1:length(d)-150));
    [~, pos_peak_locations] = max(d(1:neg_peak_locations));
    
    % One of the APA peak is between the both crossover calculated above.
    % The second one is the minimum value after the descent.
    [~,index_APA_X_C1_1 ]= max(comp_X_1 (pos_peak_locations: neg_peak_locations));
    index_APA_X_C1_1 = index_APA_X_C1_1+ pos_peak_locations;
    
    [~,index_APA_X_C1_2] = min(comp_X_1(neg_peak_locations:length(comp_X_1)));
    index_APA_X_C1_2 = index_APA_X_C1_2 + neg_peak_locations;
    
    % Second component.
    % We have to find the minumum of this signals. This is when the patient
    % goes forward.
    [~,index_APA_X_C2_1]= min(comp_X_2);
    
    % Stores the signals in  single variable ( in a row) to compare after
    % between patients.
    Signal_X = reshape(X_x ,[1, 3*max_length]);
    
    % Interpotate in each iteration.
    if(i-3)==1
        Signals_X (i-3,:) = Signal_X;
    else
     % Check what length is larger to interpolate afterwards.
      if (length(Signal_X)<length(Signals_X(i-4,:))) % If the new signal is smaller, we interpolate this.
         Signal_X =  interp1([1:length(Signal_X)],Signal_X,[1:length(Signals_X(i-4,:))]);
         Signals_X (i-3,:) = Signal_X;

      else % If the new signal is higher, we interpolate the rest of the signals.
          [r,c] = size(Signals_X);
          Signals_aux = Signals_X;
          Signals_X = zeros (1,length(Signal_X));
          
          for n= 1:r
              Signals_X (n,:) = interp1([1:c],Signals_aux (n,:),[1:length(Signal_X)]);
          end
          Signals_X (i-3,:) = Signal_X;
      end
    end
    
% -------------------------------------------------------------------------
% 4.2) PCA in Medio-Lateral Direction.
% -------------------------------------------------------------------------
    % Interpolation of the signals.
    max_length = max([length(ML_COP_mean),length(a_trunk_data_Y_mean),...
                length(g_trunk_data_Y_mean)]);

    ML_COP_mean = interp1([1:length(ML_COP_mean)],ML_COP_mean,[1:max_length]);
    a_trunk_data_Y_mean = interp1([1:length(a_trunk_data_Y_mean)],...
                            a_trunk_data_Y_mean,[1:max_length]);
    g_trunk_data_Y_mean = interp1([1:length(g_trunk_data_Y_mean)],...
                            g_trunk_data_Y_mean,[1:max_length]);
    % Center the data.
    ML_COP_mean_c = ML_COP_mean - mean(find(ML_COP_mean~=NaN));
    a_trunk_data_Y_mean_c = a_trunk_data_Y_mean - mean(find(a_trunk_data_Y_mean~=NaN));
    g_trunk_data_Y_mean_c = g_trunk_data_Y_mean - mean(find(g_trunk_data_Y_mean~=NaN));

    % Apply PCA.
    X_y = [ML_COP_mean_c; a_trunk_data_Y_mean_c; g_trunk_data_Y_mean_c]';

    [COEFF_y,SCORE_y,latent,tsquare] = princomp(X_y,'econ');
    
    % The first componets of the SCORE are the most variability, so we will
    % use them to extract features with more accuracity.
    comp_Y_1 = SCORE_y(:,1);
    comp_Y_2 = SCORE_y(:,2);
    
    % Feature extration of the componentes in ML direction.
    % First component.
    % Detect when there is a strong change of level in the signal.
     d=diff(comp_Y_1);
     [~, neg_peak_locations] = min(d(1:length(d)-50));
     
    % Obtain the APA peak in ML direction. This is the maximum value before
    % this, close to the descent.
     [~, pos_peak_locations] = max(comp_Y_1(neg_peak_locations-50:neg_peak_locations));
     index_APA_Y_C1_1 = pos_peak_locations + neg_peak_locations - 50;
    
    % Second component.
    % We have to find the minumum of this signals and after that the prior
    % maximum. This is when the patient moves toward left and right.
    [~,index_APA_Y_C2_1]= min(comp_Y_2);
    
    [~,index_APA_Y_C2_2]= max(comp_Y_2(1:index_APA_Y_C2_1));
    
    % Stores the signals in  single variable ( in a row) to compare after
    % between patients.
    Signal_Y = reshape(X_y,[1, 3*max_length]);
    
    % Interpotate in each iteration.
    if(i-3)==1
        Signals_Y (i-3,:) = Signal_Y;
    else
     % Check what length is larger to interpolate afterwards.
      if (length(Signal_Y)<length(Signals_Y(i-4,:))) % If the new signal is smaller, we interpolate this.
         Signal_Y =  interp1([1:length(Signal_Y)],Signal_Y,[1:length(Signals_Y(i-4,:))]);
         Signals_Y (i-3,:) = Signal_Y;

      else % If the new signal is higher, we interpolate the rest of the signals.
          [r,c] = size(Signals_Y);
          Signals_aux = Signals_Y;
          Signals_Y = zeros (1,length(Signal_Y));
          
          for n= 1:r
              Signals_Y (n,:) = interp1([1:c],Signals_aux (n,:),[1:length(Signal_Y)]);
          end
          Signals_Y (i-3,:) = Signal_Y;
      end
    end

    % Calculate the APA Parameters for PCA signals.
    APA_X_C1_1 = abs(comp_X_1(index_APA_X_C1_1) - comp_X_1(1));
    APA_X_C1_2 = abs(comp_X_1(index_APA_X_C1_2) - comp_X_1(index_APA_X_C1_1));
    APA_X_C1_3 = abs(comp_X_1(index_APA_X_C1_2) - comp_X_1(1));

    APA_Y_C1_1 = comp_Y_1(index_APA_Y_C1_1);

    APA_X_C2_1 =  abs( comp_X_2(index_APA_X_C2_1) - comp_X_2(1));

    APA_Y_C2_1 =  abs( comp_Y_2(index_APA_Y_C2_1) - comp_Y_2(1));
    APA_Y_C2_2 =  abs( comp_Y_2(index_APA_Y_C2_1) - comp_Y_2(index_APA_Y_C2_2));
    APA_Y_C2_3  =  abs( comp_Y_2(index_APA_Y_C2_2) - comp_Y_2(1));
    
    APA_C1_duration = abs(AP_COP_complete_ts.time(index_APA_Y_C1_1) - ...
    AP_COP_complete_ts.time(initcross_COP(1)+length(comp_Y_1)-1));

    APA_C2_duration = a_trunk_complete_ts.time(index_APA_Y_C2_2) - ...
        a_trunk_complete_ts.time(index_APA_Y_C2_1);

    APA_Parameters_PCA(i-3,:) = [APA_X_C1_1, APA_X_C1_2, APA_X_C1_3, APA_Y_C1_1,...
                     APA_X_C2_1, APA_Y_C2_1, APA_Y_C2_2, APA_Y_C2_3,...
                     APA_C1_duration, APA_C2_duration];
 end
end


%--------------------------------------------------------------------------
% 5) Correlations.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 5.1) Correlations in the original signals.
%--------------------------------------------------------------------------
% Acc and COP.
[corr_AP1, ~] = corr(APA_Parameters(:,1),APA_Parameters(:,5)); 
[corr_AP2, ~] = corr(APA_Parameters(:,2),APA_Parameters(:,5));
[corr_AP3, ~] = corr(APA_Parameters(:,3),APA_Parameters(:,5));

[corr_ML1,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,6));
[corr_ML2,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,7));
[corr_ML3,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,8));

% Gyro and COP.
[corr_APg1, ~] = corr(APA_Parameters(:,1),APA_Parameters(:,9)); 
[corr_APg2, ~] = corr(APA_Parameters(:,2),APA_Parameters(:,9));
[corr_APg3, ~] = corr(APA_Parameters(:,3),APA_Parameters(:,9));

[corr_MLg1,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,10));
[corr_MLg2,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,11));
[corr_MLg3,~] =  corr(APA_Parameters(:,4),APA_Parameters(:,12));

% Duration.
[corr_dur1,~] =  corr(APA_Parameters(:,13),APA_Parameters(:,14));
[corr_dur2,~] =  corr(APA_Parameters(:,14),APA_Parameters(:,15));
[corr_dur3,~] =  corr(APA_Parameters(:,13),APA_Parameters(:,15));

% Lump together the diferents features calculated above.
% In Antero-Posterior direction.
% 1-Corr between positive peak in AP-COP and the negative peak in the Acc.
% 2-Corr between the peaks distance in AP-COP and the negative peak in the Acc.
% 3-Corr between negative peak in AP-COP and the negative peak in the Acc.
% 4-Corr between positive peak in AP-COP and the negative peak in the Gyro.
% 5-Corr between the peaks distance in AP-COP and the negative peak in the Gyro.
% 6-Corr between negative peak in AP-COP and the negative peak in the Gyro.
correlationsAP = [corr_AP1, corr_AP2, corr_AP3, corr_APg1, corr_APg2,...
                    corr_APg3];
                
% In Medio-Lateral direction.
% 1-Corr between positive peak in ML-COP and the positive peak in the Acc.
% 2-Corr between positive peak in ML-COP and the peaks distance in the Acc.
% 3-Corr between positive peak in ML-COP and the negative peak in the Acc.
% 4-Corr between positive peak in ML-COP and the positive peak in the Gyro.
% 5-Corr between positive peak in ML-COP and the peaks distance in the Gyro.
% 6-Corr between positive peak in ML-COP and the negative peak in the Gyro.
correlationsML = [corr_ML1, corr_ML2, corr_ML3, corr_MLg1, corr_MLg2,...
                    corr_MLg3];
% APA Duration.
% 1-Corr between AP-COP duration and Acc-duration.
% 2-Corr between Acc-duration and Gyro duration.
% 3-Corr between AP-COP duration and Gyro duration.
correlationsDuration = [corr_dur1,corr_dur2,corr_dur3];
                
% -------------------------- Correlations Plots----------------------------
if strcmpi(showPlotsCorr,'yes')
    
figure()              
x_axes={'Acc-COP_AP1','Acc-COP_AP2','Acc-COP_AP3',...
    'Gyro-COP_AP1','Gyro-COP_AP2','Gyro-COP_AP3'};

bar (correlationsAP);
set(gca,'XtickL',x_axes)
title('Differents Correlations between measures of APAs in AP direction'); 

figure()              
x_axes={'Acc-COP_ML1','Acc-COP_ML2','Acc-COP_ML3',...
    'Gyro-COP_ML1','Gyro-COP_ML2','Gyro-COP_ML3'};

bar (correlationsML);
set(gca,'XtickL',x_axes)
title('Differents Correlations between measures of APAs in ML direction');

figure()              
x_axes={'Dur_COP-Acc','Dur_Acc-Gyro','Dur_COP-Gyro'};

bar (correlationsDuration);
set(gca,'XtickL',x_axes)
title('Differents Correlations between measures of APAs duration');

end
fprintf(' Feature extraction completed !!! \n');

%--------------------------------------------------------------------------
% 5.2) Correlations in the projections of the signals after applying PCA.
%--------------------------------------------------------------------------

% First component.
[corr_X1, ~] = corr(APA_Parameters_PCA(:,1),APA_Parameters_PCA(:,5)); 
[corr_X2, ~] = corr(APA_Parameters_PCA(:,2),APA_Parameters_PCA(:,5));
[corr_X3, ~] = corr(APA_Parameters_PCA(:,3),APA_Parameters_PCA(:,5));

[corr_Y1,~] =  corr(APA_Parameters_PCA(:,4),APA_Parameters_PCA(:,6));
[corr_Y2,~] =  corr(APA_Parameters_PCA(:,4),APA_Parameters_PCA(:,7));
[corr_Y3,~] =  corr(APA_Parameters_PCA(:,4),APA_Parameters_PCA(:,8));

[corr_duration,~] = corr(APA_Parameters_PCA(:,9),APA_Parameters_PCA(:,10));

% Lump together the diferents features calculated above.
% In Antero-Posterior direction.
% 1-Corr between positive peak in comp1 and the negative peak in comp2.
% 2-Corr between the peaks distance in comp1 and the negative peak in the
% comp2
% 3-Corr between negative peak in comp1 and the negative peak in comp2.
% In Medio-Lateral direction.
% 4-Corr between positive peak in comp1 and the positive peak in comp2.
% 5-Corr between positive peak in comp1 and the peaks distance in comp2.
% 6-Corr between positive peak in comp1 and the negative peak in comp2.
% 7- Corr between comp1-APA duration and comp2-APA duration.
correlationsPCA = [corr_X1, corr_X2, corr_X3,corr_Y1, corr_Y2, corr_Y3,...
    corr_duration];

% Plots
if strcmpi(showPlotsCorr,'yes')   
figure()              
x_axes={'Comp_X1','Comp_X2','Comp_X3',...
    'Comp_Y1','Comp_Y2','Comp_Y3', 'Durat'};

bar (correlationsPCA);
set(gca,'XtickL',x_axes)
title('Differents Correlations between components after applying PCA'); 
end
%--------------------------------------------------------------------------
% 6) PCA between patients.
%--------------------------------------------------------------------------
% Apply PCA.
Signals_X = Signals_X - mean(find(Signals_X~=NaN));
Signals_Y = Signals_Y - mean(find(Signals_Y~=NaN));

[COEFF_X,SCORE_X,latent,tsquare] = princomp(Signals_X','econ');
[COEFF_Y,SCORE_Y,latent,tsquare] = princomp(Signals_Y','econ');

figure()
biplot(COEFF_X(:,1:2),'Scores',SCORE_X(:,1:2));
title('PCA between patients in Antero-Posterior direction')

figure()
biplot(COEFF_Y(:,1:2),'Scores',SCORE_Y(:,1:2));
title('PCA between patients in Medio-Lateral direction')
