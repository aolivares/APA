function [ COP_AP_total, COP_ML_total] = Calculation_COP( folder)
% CALCULATION_COP calculate the center of pressure antero-posterior and 
% medio-lateral of all of data in a folder.
%
% - Input:
%    |_ 'folder': folder where there are all of data to calculate the COP.
   
% - Output:
%    |_ ' COP_AP_total': Antero-Posterior Center of Pressure  of all data
%    in the folder.
%    |_ ' COP_ML_total': Medio-Lateral Center of Pressure  of all data
%    in the folder.
%
% -------------------------------------------------------------------------

% Read all datafiles
fnames = dir([folder,'*.txt']);
numfids = length(fnames);

for n = 1:numfids
  namefile = fnames(n).name;
  
  % Read formatted data from text file. As always "fid" is a file identifier 
  % that we obtain with "fopen".
  fid = fopen(fullfile(folder,namefile));
  

    % Data is a matrix of 19 columns with each parameter. 'repmat' is  function 
    % to indicate that the %n (double) conversion specifier should appear 19 
    % times. This technique is useful when a format repeats many times. 
    data = textscan(fid,repmat('%n',[1,19]),'CollectOutput',1);
    data = data{1,1};
    
    % The format of the data is (by columns):
    % 1 : time
    % 2-9 : force of echa sensor (8 sensors) of the left foot (vertical ground
    % reaction).
    % 10-17 : force of each sensors of the right foot.
    % 18 : sum left force.
    % 19 : sum right force.

    % Obtain the variables.
    sample_f =100;
    time = data(:,1);
    num_samples = length(time);
    force_sensors_left = data(:,2:9);
    force_sensors_right = data(:,10:17);
    force_left = data (:,18);
    force_right = data (:,19);
    force_total = sum([force_left, force_right]');

    %----------------------------------------------------------------------
    % 2) Calculate the center of pressure.
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % 2.1) Calculate ML-COP.
    %----------------------------------------------------------------------

    % Obtain the force in the ML-direction. (To undertand this, seet he
    % coordenated of each sensors).
    % we calculate the sum of the force of the aligned sensors in the
    % medio-laretal direction (same coordenates in the X axis).

    force_ML(:,1) = force_sensors_left(:,2) + force_sensors_left(:,4) + ...
        force_sensors_left(:,6);

    force_ML(:,2) = force_sensors_left(:,1) + force_sensors_left(:,8);

    force_ML(:,3) = force_sensors_left(:,7) + force_sensors_left(:,5) + ...
        force_sensors_left(:,3);

    force_ML(:,4) = force_sensors_right(:,3) + force_sensors_right(:,5) + ...
        force_sensors_right(:,7);

    force_ML(:,5) = force_sensors_right(:,1) + force_sensors_right(:,8);

    force_ML(:,6) = force_sensors_right(:,6) + force_sensors_right(:,4) + ...
        force_sensors_right(:,2);

    coord_X =[-700, -500, -300, 300, 500, 700];

    % Calculate the center of pressure.
    for i=1:num_samples
        COP_ML(i) = sum(force_ML(i,:).*coord_X)/sum(force_ML(i,:));
    end

    %----------------------------------------------------------------------
    % 2.2) Calculate AP-COP.
    %----------------------------------------------------------------------

    % Obtain the force in the AP-direction. (To undertand this, seet he
    % coordenated of each sensors).
    % We calculate the sum of the force of the aligned sensors in the
    % antero-posterior direction (same coordenates in the Y axis).

    force_AP(:,1) = force_sensors_left(:,1) + force_sensors_right(:,1);

    force_AP(:,2) = force_sensors_left(:,2) + force_sensors_left(:,3) + ...
                        force_sensors_right(:,2) + force_sensors_right(:,3);

    force_AP(:,3) = force_sensors_left(:,4) + force_sensors_left(:,5) + ...
                        force_sensors_right(:,4) + force_sensors_right(:,5);

    force_AP(:,4) =  force_sensors_left(:,6) + force_sensors_left(:,7) + ...
                        force_sensors_right(:,6) + force_sensors_right(:,7);

    force_AP(:,5) = force_sensors_left(:,8) + force_sensors_right(:,8);

    coord_Y =[-800, -400, 0, 400, 800];

    % Calculate the center of pressure.
    for i=1:num_samples
        COP_AP(i) = sum(force_AP(i,:).*coord_Y)/sum(force_AP(i,:));
    end

    %----------------------------------------------------------------------
    % 3) Obtain a single signal to characterise the person.
    %----------------------------------------------------------------------
    % We have to detect the differents steps that the person does during the
    % experiment and calculate the mean. Thus, we will a have a single signal
    % of ML and AP COP that define the movement of the person.

    %----------------------------------------------------------------------
    % 3.1) Determine the bounds of each step.
    %----------------------------------------------------------------------
    % We are going to delimit the intervals detecting the point when the person
    % is lifting one of the foot and going down the other. We can detect this
    % choosing one of the signal, for example, right foot and determining the
    % samples where the force changes from the low value to high value, i.e
    % when there are a positive peak in the derivate of the signal.
    d = diff(force_right);
    [peak_values, peak_locations] = findpeaks(d, 'minpeakheight', 20, ...
        'MinPeakDistance',90);
    
    %----------------------------------------------------------------------
    % 3.2) Align signals and carry out the average.
    %----------------------------------------------------------------------
    COP_AP_mean= aligned_signals( COP_AP(peak_locations(1):peak_locations(2)),...
        COP_AP(peak_locations(2):peak_locations(3)) );

    COP_ML_mean= aligned_signals( COP_ML(peak_locations(1):peak_locations(2)),...
        COP_ML(peak_locations(2):peak_locations(3)) );

    for i=3:length(peak_locations)-1
        COP_AP_mean = aligned_signals(COP_AP(peak_locations(i):peak_locations(i+1)),...
            COP_AP_mean);

        COP_ML_mean = aligned_signals(COP_ML(peak_locations(i):peak_locations(i+1)),...
            COP_ML_mean);
    end

    
    % Interpotate in each iteration.
    if n==1
            
    COP_AP_total (n,:) = COP_AP_mean;
    COP_ML_total (n,:) = COP_ML_mean;
    
    else
        
     % Check what length is larger to interpolate afterwards.
     % We check only one of the signals because they have the same length.
      if (length(COP_AP_mean)<length(COP_AP_total(n-1,:))) % If the new signal is smaller, we interpolate this.
         COP_AP_mean =  interp1([1:length(COP_AP_mean)],COP_AP_mean,...
                        [1:length(COP_AP_total(n-1,:))]);
         COP_ML_mean =  interp1([1:length(COP_ML_mean)],COP_ML_mean,...
                        [1:length(COP_ML_total(n-1,:))]);
         
         COP_AP_total (n,:) = COP_AP_mean;
         COP_ML_total (n,:) = COP_ML_mean;

      else % If the new signal is higher, we interpolate the rest of the signals.
          [r_AP,c_AP] = size(COP_AP_total);
          [r_ML,c_ML] = size(COP_ML_total);
          
          COP_AP_aux = COP_AP_total;
          COP_ML_aux = COP_ML_total;
          
          COP_AP_total = zeros (1,length(COP_AP_mean));
          COP_ML_total = zeros (1,length(COP_ML_mean));
          
          for k= 1:r_AP
              COP_AP_total  (k,:) = interp1([1:c_AP],COP_AP_aux (k,:),...
                                    [1:length(COP_AP_mean)]);
          end
          for k= 1:r_ML
              COP_ML_total  (k,:) = interp1([1:c_ML],COP_ML_aux (k,:),...
                                    [1:length(COP_ML_mean)]);
          end
          
          COP_AP_total (n,:) = COP_AP_mean;
          COP_ML_total (n,:) = COP_ML_mean;
      end
    end
end

end

