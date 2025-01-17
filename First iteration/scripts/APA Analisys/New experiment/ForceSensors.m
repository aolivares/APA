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
folder_Co = 'database/Control/';
[COP_AP_Co, COP_ML_Co] = Calculation_COP (folder_Co);

folder_Pt = 'database/Patients/';
[COP_AP_Pt, COP_ML_Pt] = Calculation_COP (folder_Pt);

% Obtain the lenghts of the signals to interpolate because we need all of
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

% Match both signals that represents the person and group all together.
X_Co = [COP_AP_Co_int COP_ML_Co_int]';
X_Pt = [COP_AP_Pt COP_ML_Pt]';

% Apply PCA for all data.
% Group all data together.
X = [X_Co X_Pt];

[COEFF,SCORE,latent,tsquare, explained] = princomp(X,'econ');
figure()
biplot(COEFF(:,1:2),'Scores',SCORE(:,1:2));

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
c = zeros(1,length(X_Co(1,:)));
p = ones(1,length(X_Pt(1,:)));
labels=[c p];

%--------------------------------------------------------------------------
% 4) SVM algorithm for clasification.
%--------------------------------------------------------------------------
% Apply PLS for differents components (latents variables) and SVM for each
% case.
for i=1:ncomp

    % Apply PLS.
    XS = PLS_feature_extraction2(labels',X',i);
    
    % Apply SVM for classification.
    % Obtain the parameters to apply the algorithm.
    labels = logical (labels);
    [acc, sen, spe] = SVM_LOU(XS,labels);
    acc_pls (i,:) = acc;
    sen_pls (i,:) = sen;
    spe_pls (i,:) = spe;
    
end

% SVM for row data.
[acc_row, sen_row, spe_row] = SVM_LOU(X',labels);

% SVM for PCA
for i= 1:4
    [acc, sen, spe] = SVM_LOU(SCORE(1:i,:)',labels);
    acc_pca(i,:) = acc;
    sen_pca (i,: ) = sen;
    spe_pca (i,:) = spe;
end

% Figures
% Make a scree plot of the percent variability explained by each principal
% component.
figure()
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

% Stadistics parameters for each number of components for PLS algorithm.
figure()
subplot(3,1,1)
plot(acc_pls)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Accuracy(%)')
title('Accuracity with PLS algorithm')

subplot(3,1,2)
plot(sen_pls)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Sensitivity')
title('Sensitivity with PLS algorithm')

subplot(3,1,3)
plot(spe_pls)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Specificity')
title('Specificity with PLS algorithm')

% Stadistics parameters for each number of components for PCA algorithm.
figure()
subplot(3,1,1)
plot(acc_pca)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Accuracy(%)')
title('Accuracity with PCA algorithm')

subplot(3,1,2)
plot(sen_pca)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Sensitivity')
title('Sensitivity with PCA algorithm')

subplot(3,1,3)
plot(spe_pca)
legend('linear','quadratic','polynomial','rbf')
xlabel('Number of components')
ylabel('Specificity')
title('Specificity with PCA algorithm')

% Plot 3D of the classification
NORMAL = 1:18;
DTA = 19:45;

% Classifiction with PLS (kern = lineal).
figure()
y= draw_dec_surf_svm(XS, NORMAL, DTA, 'linear');

% Classification with PLS (kern = rbf)
figure()
y= draw_dec_surf_svm(XS, NORMAL, DTA, 'polynomial');

% Classification with PCA (pral comp=4, kern = lineal)
figure()
y= draw_dec_surf_svm(SCORE(1:3,:)', NORMAL, DTA, 'rbf');

% Apply ROC.
% ROC curve for linear kernel and seven components.
 XS = PLS_feature_extraction2(labels',X',7);
 step = 500;
 
figure();
[X_PLS7,Y_PLS7,auc_PLS7]=roc(XS,labels,step);

% ROC curve for linear kernel in PCA data.
figure()
[X_PCA,Y_PCA,auc_PCA]=roc(SCORE(1:4,:)',labels,step);

% Comparative representation of ROC.
figure()
plot(X_PLS7,Y_PLS7,'b');
hold on;
plot(X_PCA,Y_PCA,'g');

% -----------------------------------------------------------------------
% Plot PCA, PLS and original mean signals.
X = X';
Y = labels';
ncomp = 7;
P= size(X,1);
XS= zeros(P,ncomp);
meanX = mean(X,1);
X0 = bsxfun(@minus, X, meanX);

for p=1:P
%     disp(['PLS feature extraction for patient ' num2str(p) ]);
    train = true(P,1);
    train(p)= false;
    [XL,yL,XSp,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(train,:),Y(train),ncomp); 
    W= stats.W;
    XS(p,:)= X0(p,:)*W;
end

X_pls = XS*XL';
X_pls_mean = sum(X_pls);

X_pca_mean = sum(SCORE(:,4),2);

X_mean = mean(X);

% Replace the NaN values per 0s for PCA data.
for i =1: length(X_pca_mean)
       if( isnan(X_pca_mean (i,1))==1)
            X_pca_mean (i,1)= 0;
       end
end

figure;
subplot(3,1,1)
plot(X_mean)

subplot (3,1,2)
plot(X_pca_mean)

subplot(3,1,3)
plot(X_pls_mean)

