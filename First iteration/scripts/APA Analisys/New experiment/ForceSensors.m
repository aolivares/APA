% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------- Classification force data -------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      
%                 - Veronica  Torres              
%                 - Dr. Eng. Alberto Olivares 
%                 - Dr. Eng. Juan Manuel Gorriz 
%
% * Last modification:22/06/2015

% INFORMATION: 
% 
% * 1) Trainning
% 
% * 2) Testing.

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Trainning.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Obtain the data.
%--------------------------------------------------------------------------

% Delete workspace.
clear all, close all, clc;

% Calculate a single signal per patient.
folder_Co = 'database/trainning/Control/';
[COP_AP_Co, COP_ML_Co] = Calculation_COP (folder_Co);

folder_Pt = 'database/trainning/Patients/';
[COP_AP_Pt, COP_ML_Pt] = Calculation_COP (folder_Pt);

L_AP = length(COP_AP_Pt);
L_ML = length(COP_ML_Pt);
% Interpolate the signals from control subject for these have the same
% length.
for i =1: length(COP_AP_Co(:,1))
    COP_AP_Co_int(i,:) = interp1([1:length(COP_AP_Co(1,:))],COP_AP_Co(i,:),[1:length(COP_AP_Pt(1,:))]);
    
    COP_ML_Co_int(i,:) = interp1([1:length(COP_ML_Co(1,:))],COP_ML_Co(i,:),[1:length(COP_ML_Pt(1,:))]);
    
end



% -------------------------------------PCA---------------------------------
% Apply PCA for control subjects.
X_Co = [COP_AP_Co_int COP_ML_Co_int]';

% Center the data
for i =1:length(X_Co(1,:))
    points_mean = find(abs(X_Co(:,i))>0);
    X_Co (:,i) = X_Co(:,i) - mean(X_Co(points_mean,i));

end

[COEFF_Co,SCORE_Co,latent,tsquare] = princomp(X_Co,'econ');
Comp_Co = SCORE_Co(:,1);


% Apply PCA for Patients.
X_Pt = [COP_AP_Pt COP_ML_Pt]';

% Center the data
for i =1:length(X_Pt(1,:))
    points_mean = find(abs(X_Pt(:,i))>0);
    X_Pt (:,i) = X_Pt(:,i) - mean(X_Pt(points_mean,i));
end
[COEFF_Pt,SCORE_Pt,latent,tsquare] = princomp(X_Pt,'econ');
Comp_Pt = SCORE_Pt(:,1);

% Plot
figure; plot(Comp_Co); hold on; plot(Comp_Pt,'g');

% Apply PCA for all data
X = [X_Co X_Pt];

[COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
Comp_1= SCORE(:,1);
Comp_2 = SCORE (:,2);


% Determine the threshold using the eclude distance between the components
% calculated.
% For control subjects:
distance_1 = abs(Comp_1 - Comp_Co);
distance_1 = mean(distance_1(find(distance_1>0)));
distance_2 = abs(Comp_2 - Comp_Co);
distance_2 = mean(distance_2(find(distance_2>0)));

threshold_Co = min([distance_1,distance_2]);

% For patients
distance_1 = abs(Comp_1 - Comp_Pt);
distance_1 = mean(distance_1(find(distance_1>0)));
distance_2 = abs(Comp_2 - Comp_Pt);
distance_2 = mean(distance_2(find(distance_2>0)));

threshold_Pt = min([distance_1,distance_2]);

% Clear variables.
clearvars -except threshold_Pt threshold_Co Comp_Pt Comp_Co X_Co X_Pt L_ML ...
L_AP

%--------------------------------------------------------------------------
% Testing
%--------------------------------------------------------------------------
% Calculate a single signal per patient.
folder_Co = 'database/testing/Control/';
[COP_AP_Co, COP_ML_Co] = Calculation_COP (folder_Co);

folder_Pt = 'database/testing/Patients/';
[COP_AP_Pt, COP_ML_Pt] = Calculation_COP (folder_Pt);

% Interpolate the control subjects.
for i =1: length(COP_AP_Co(:,1))
    COP_AP_Co_int(i,:) = interp1([1:length(COP_AP_Co(1,:))],COP_AP_Co(i,:),[1:L_AP]);
    
    COP_ML_Co_int(i,:) = interp1([1:length(COP_ML_Co(1,:))],COP_ML_Co(i,:),[1:L_ML]);
    
end

% Interpolate the patients.
for i =1: length(COP_AP_Pt(:,1))
    COP_AP_Pt_int(i,:) = interp1([1:length(COP_AP_Pt(1,:))],COP_AP_Pt(i,:),[1:L_AP]);
    
    COP_ML_Pt_int(i,:) = interp1([1:length(COP_ML_Pt(1,:))],COP_ML_Pt(i,:),[1:L_ML]);
    
end

COP_Co = [COP_AP_Co_int COP_ML_Co_int];
COP_Pt = [COP_AP_Pt_int COP_ML_Pt_int];

% For control subject
for i = 1:length(COP_Co(:,1))
    
    % Center the data.
    COP_CO = COP_Co(i,:);
    points_mean = find(abs(COP_CO)>0);
    COP_CO = COP_CO - mean(COP_CO(points_mean));
    
    % Apply PCA.
    X = [X_Co COP_CO'];
    [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
    Comp_1  = SCORE(:,1);
    
    % Calcule the distance.
    distance = abs(Comp_1 - Comp_Co);
    distance_Co_Co (i) = mean(distance(find(distance>0)));
    
    % Apply PCA.
     X = [X_Pt COP_CO'];
    [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
    Comp_2 = SCORE(:,1);
    
    % Calcule the distance.
    distance = abs(Comp_2 - Comp_Pt);
    distance_Co_Pt (i) = mean(distance(find(distance>0)));
    
end

% For patients
for i = 1:length(COP_Pt(:,1))

    COP_PT = COP_Pt(i,:);
    
    % Center the data.
    points_mean = find(abs(COP_PT)>0);
    COP_PT = COP_PT - mean(COP_PT(points_mean));
    
    % Apply PCA.
    X = [X_Pt COP_PT'];
    [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
    Comp_1 = SCORE(:,1);
    
    % Calcule the distance.
    distance = abs(Comp_1 - Comp_Pt);
    distance_Pt_Pt (i) = mean(distance(find(distance>0)));
    
    % Apply PCA.
     X = [X_Co COP_PT'];
    [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
    Comp_2 = SCORE(:,1);
    
    % Calcule the distance.
    distance = abs(Comp_2 - Comp_Co);
    distance_Pt_Co (i) = mean(distance(find(distance>0)));
    
end

find(distance_Co_Co<distance_Co_Pt)
find(distance_Pt_Pt<distance_Pt_Co)






% % Obtain the patter for each people.
% COP_Co = [COP_AP_Co COP_AP_Co];
% COP_Pt = [COP_AP_Pt COP_ML_Pt];
% 
% % Calcule the distance between every new people
% % Interpolate the signals from control subject for these have the same
% % length.
% for i =1: length(COP_Co(:,1))
%     COP_Co_int(i,:) = interp1([1:length(COP_Co(1,:))],COP_Co(i,:),[1:length(Comp_Co)]);
%     
%     % Calculate the distance.
%     distance = abs(COP_Co_int(i,:) - Comp_Co');
%     distance_Co_Co (i) = mean(distance(find(distance>0)));
%     
%     distance = abs(COP_Co_int(i,:) - Comp_Pt');
%     distance_Co_Pt (i) = mean(distance(find(distance>0)));
% end
% 
% for i =1: length(COP_Pt(:,1))
%     COP_Pt_int(i,:) = interp1([1:length(COP_Pt(1,:))],COP_Pt(i,:),[1:length(Comp_Pt)]);
%     
%     % Calculate the distance.
%     distance = abs(COP_Pt_int(i,:) - Comp_Co');
%     distance_Pt_Co (i) = mean(distance(find(distance>0)));
%     
%     distance = abs(COP_Pt_int(i,:) - Comp_Pt');
%     distance_Pt_Pt (i) = mean(distance(find(distance>0)));
% end
% 
% % Identify the class of each person.
% Co_ok = length(find(distance_Co_Co<distance_Co_Pt))/length(distance_Co_Co)*100;
% Pt_ok = length(find(distance_Pt_Pt<distance_Pt_Co))/length(distance_Pt_Co)*100;
