% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------- Classification force data -------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data. 
%
% * Authors:      
%                 - Veronica  Torres (1) 
%                   |_ vts24@correo.ugr.es
%                 - Dr. Eng. Alberto Olivares (2)
%                   |_ aolivares@ugr.es
%                 - Dr. Eng. Juan Manuel Gorriz (3)
%                   |_ gorriz@ugr.es
%
% * Affiliation: (1) Master in Telecommunication Engineering, University of 
%                    Granada, Granada, Spain, (student).
%                (2) Signal Processing and Biomedical Applications group,
%                    Department of Signal Theory, Telematics and
%                    Communications, University of Granada, Granada, Spain
%                (3) Signal Processing and Biomedical Applications group,
%                    Department of Signal Theory, Telematics and
%                    Communications, University of Granada, Granada, Spain.
%
% * Last modification: 25/06/2015

% INFORMATION: This file contains the routine to characterise and clasify
% the gait wile people (patients and control subject) walk normally.
% 
% * 1) Obtain data.
% 
% * 2) Apply PCA to extract features.
% 
% * 3) PLS feature extraction.
%
% * 4) SVM algorithm for clasification.
%
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
% Delete workspace.
clear all, close all, clc;

%--------------------------------------------------------------------------
% 1) Obtain the data.
%--------------------------------------------------------------------------
% Calculate the Center of Pressure of all data that corresponds to one step
% for each people. Thus, we have two signals per people (AP-COP and
% ML-COP).

% Calculate a single signal per patient.
folder_Co = 'database/trainning/Control/';
[COP_AP_Co, COP_ML_Co] = Calculation_COP (folder_Co);

folder_Pt = 'database/trainning/Patients/';
[COP_AP_Pt, COP_ML_Pt] = Calculation_COP (folder_Pt);

% Obtain the lenghts os the signals to interpolate because we need all of
% them have the same dimention.
L_AP = length(COP_AP_Pt);
L_ML = length(COP_ML_Pt);

% Interpolate the signals from control subject for these have the same
% length.
for i =1: length(COP_AP_Co(:,1))
    COP_AP_Co_int(i,:) = interp1([1:length(COP_AP_Co(1,:))],COP_AP_Co(i,:),...
                        [1:L_AP]);
    
    COP_ML_Co_int(i,:) = interp1([1:length(COP_ML_Co(1,:))],COP_ML_Co(i,:),...
                        [1:L_ML]);
    
end


%--------------------------------------------------------------------------
% 2) Apply PCA to extract features.
%--------------------------------------------------------------------------

% Apply PCA for control subjects.
% Match both signals that represents the person and group all together.
X_Co = [COP_AP_Co_int COP_ML_Co_int]';

% % Center the data
% for i =1:length(X_Co(1,:))
%     points_mean = find(abs(X_Co(:,i))>0);
%     X_Co (:,i) = X_Co(:,i) - mean(X_Co(points_mean,i));
% 
% end
% 
% % Apply PCA.
% [COEFF_Co,SCORE_Co,latent,tsquare] = princomp(X_Co,'econ');
% Comp_Co = SCORE_Co(:,1);


% Apply PCA for Patients.
X_Pt = [COP_AP_Pt COP_ML_Pt]';

% % Center the data
% for i =1:length(X_Pt(1,:))
%     points_mean = find(abs(X_Pt(:,i))>0);
%     X_Pt (:,i) = X_Pt(:,i) - mean(X_Pt(points_mean,i));
% end
% 
% % Apply PCA.
% [COEFF_Pt,SCORE_Pt,latent,tsquare] = princomp(X_Pt,'econ');
% Comp_Pt = SCORE_Pt(:,1);

% Plot
% figure; plot(Comp_Co); hold on; plot(Comp_Pt,'g');

% Apply PCA for all data.
% Group all data together.
X = [X_Co X_Pt];

[COEFF,SCORE,latent,tsquare] = princomp(X,'econ');

%--------------------------------------------------------------------------
% 3) PLS feature extraction.
%--------------------------------------------------------------------------

% Replace the NaN values per 0s for PLS can work properly.
for i =1: length(X(1,:))
    for j= 1:length(X(:,1))
       if( isnan(X (j,i))==1)
            X (j,i)= 0;
       end
    end
end

% Establish the parameters for PLS regresion.
% ncomp : latents variables.
% labels : 0--> control subjects and 1--> patients.
% X : all data.
% XS : PLs components that are linear combination of the variables in X.
ncomp=31;
c = zeros(1,14);
p = ones(1,19);
labels=[c p]';
XS = PLS_feature_extraction2(labels,X',ncomp);

%--------------------------------------------------------------------------
% 4) SVM algorithm for clasification.
%--------------------------------------------------------------------------

% Obtain the parameters to apply the algorithm.
labels = logical (labels);
data = XS;
P = size(data,1);

% Number of time that the experiment is carried out.(Leave-M-Out).
N = P;  

% Number of patients and control subjects.
NPat = sum(labels==1);
NCont = sum(labels~=1);

% Define the kinds of kernel functions.
kern = {'linear','quadratic','polynomial','rbf'};

% Obtain the stadistics for each kernel
for k = 1:length(kern)
    kernel = char(kern(k));
    Errors = ones(N,1);
    
    for n = 1:N
        
        % Define the set of trainning and testing for each iteration.
        % (Leave-M-Out)In this case one them is out.
        train_set = true(1,P); train_set(n) = false; test_set =~train_set;
        train_labels = labels(train_set); test_labels = labels(test_set);
        
        tr_data = data(train_set,:);
        te_data = data(test_set,:);
        
        % Appy SVM.
        svmStruct = svmtrain(tr_data, train_labels, 'Kernel_Function', ...
                    kernel,'Showplot', false);
        class = svmclassify(svmStruct, te_data)';
        
        
        % Leave-M-Out.
        class = cast(class,'double');
        Errors(n) = length(find(class ~= test_labels'));
        
    end
    
    
    % Leave-M-Out. Calculate the stadistics.
    Total_Errors = sum(Errors,1);
    Total_test = N;
    
    % Percentage (%) of data clasified properly.
    acc(k) = (Total_test - Total_Errors)./Total_test * 100;
    
    POS = labels==1;
    NEG = labels~=1;
    FP = sum(Errors(POS,:),1); TP = NPat-FP;
    FN = sum(Errors(NEG,:),1); TN = NCont-FN;
    
    % sen : Percentage of patients clasidied properly over the total of errors.
    % pen : Percentage of CS clasidied properly over the total of errors.
    sen(k) = (TP./(TP+FN));
    spe(k) = (TN./(TN+FP));
    
end

% oldorder = get(gcf,'DefaultAxesColorOrder');
% set(gcf,'DefaultAxesColorOrder',jet(33));
% plot3(repmat(1:33,331,1)',1:331,SCORE');
% set(gcf,'DefaultAxesColorOrder',oldorder);
% xlabel('Wavelength Index'); ylabel('Octane'); axis('tight');
% grid on


% % Determine the threshold using the eclude distance between the components
% % calculated.
% % For control subjects:
% distance_1 = abs(Comp_1 - Comp_Co);
% distance_1 = mean(distance_1(find(distance_1>0)));
% distance_2 = abs(Comp_2 - Comp_Co);
% distance_2 = mean(distance_2(find(distance_2>0)));
% 
% threshold_Co = min([distance_1,distance_2]);
% 
% % For patients
% distance_1 = abs(Comp_1 - Comp_Pt);
% distance_1 = mean(distance_1(find(distance_1>0)));
% distance_2 = abs(Comp_2 - Comp_Pt);
% distance_2 = mean(distance_2(find(distance_2>0)));
% 
% threshold_Pt = min([distance_1,distance_2]);

% % Clear variables.
% clearvars -except threshold_Pt threshold_Co Comp_Pt Comp_Co X_Co X_Pt L_ML ...
% L_AP
% 
% %--------------------------------------------------------------------------
% % Testing
% %--------------------------------------------------------------------------
% % Calculate a single signal per patient.
% folder_Co = 'database/testing/Control/';
% [COP_AP_Co, COP_ML_Co] = Calculation_COP (folder_Co);
% 
% folder_Pt = 'database/testing/Patients/';
% [COP_AP_Pt, COP_ML_Pt] = Calculation_COP (folder_Pt);
% 
% % Interpolate the control subjects.
% for i =1: length(COP_AP_Co(:,1))
%     COP_AP_Co_int(i,:) = interp1([1:length(COP_AP_Co(1,:))],COP_AP_Co(i,:),[1:L_AP]);
%     
%     COP_ML_Co_int(i,:) = interp1([1:length(COP_ML_Co(1,:))],COP_ML_Co(i,:),[1:L_ML]);
%     
% end
% 
% % Interpolate the patients.
% for i =1: length(COP_AP_Pt(:,1))
%     COP_AP_Pt_int(i,:) = interp1([1:length(COP_AP_Pt(1,:))],COP_AP_Pt(i,:),[1:L_AP]);
%     
%     COP_ML_Pt_int(i,:) = interp1([1:length(COP_ML_Pt(1,:))],COP_ML_Pt(i,:),[1:L_ML]);
%     
% end
% 
% COP_Co = [COP_AP_Co_int COP_ML_Co_int];
% COP_Pt = [COP_AP_Pt_int COP_ML_Pt_int];
% 
% % For control subject
% for i = 1:length(COP_Co(:,1))
%     
%     % Center the data.
%     COP_CO = COP_Co(i,:);
%     points_mean = find(abs(COP_CO)>0);
%     COP_CO = COP_CO - mean(COP_CO(points_mean));
%     
%     % Apply PCA.
%     X = [X_Co COP_CO'];
%     [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
%     Comp_1  = SCORE(:,1);
%     
%     % Calcule the distance.
%     distance = abs(Comp_1 - Comp_Co);
%     distance_Co_Co (i) = mean(distance(find(distance>0)));
%     
%     % Apply PCA.
%      X = [X_Pt COP_CO'];
%     [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
%     Comp_2 = SCORE(:,1);
%     
%     % Calcule the distance.
%     distance = abs(Comp_2 - Comp_Pt);
%     distance_Co_Pt (i) = mean(distance(find(distance>0)));
%     
% end
% 
% % For patients
% for i = 1:length(COP_Pt(:,1))
% 
%     COP_PT = COP_Pt(i,:);
%     
%     % Center the data.
%     points_mean = find(abs(COP_PT)>0);
%     COP_PT = COP_PT - mean(COP_PT(points_mean));
%     
%     % Apply PCA.
%     X = [X_Pt COP_PT'];
%     [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
%     Comp_1 = SCORE(:,1);
%     
%     % Calcule the distance.
%     distance = abs(Comp_1 - Comp_Pt);
%     distance_Pt_Pt (i) = mean(distance(find(distance>0)));
%     
%     % Apply PCA.
%      X = [X_Co COP_PT'];
%     [COEFF,SCORE,latent,tsquare] = princomp(X,'econ');
%     Comp_2 = SCORE(:,1);
%     
%     % Calcule the distance.
%     distance = abs(Comp_2 - Comp_Co);
%     distance_Pt_Co (i) = mean(distance(find(distance>0)));
%     
% end
% 
% find(distance_Co_Co<distance_Co_Pt)
% find(distance_Pt_Pt<distance_Pt_Co)






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
