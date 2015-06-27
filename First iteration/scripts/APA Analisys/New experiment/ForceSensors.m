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

% Match both signals that represents the person and group all together.
X_Co = [COP_AP_Co_int COP_ML_Co_int]';
X_Pt = [COP_AP_Pt COP_ML_Pt]';

% Apply PCA for all data.
% Group all data together.
X = [X_Co X_Pt];

[COEFF,SCORE,latent,tsquare, explained] = princomp(X,'econ');


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
labels=[c p]';

%--------------------------------------------------------------------------
% 4) SVM algorithm for clasification.
%--------------------------------------------------------------------------
% Apply PLS for differents components (latents variables) and SVM for each
% case.
for i=1:ncomp

    % Apply PLS.
    XS = PLS_feature_extraction2(labels,X',i);
    
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
kernel_function='linear';
% y= draw_dec_surf_svm(XS, NORMAL, DTA, 'linear');
% Figura con el hiperplano de decisión.
% Figura con el hiperplano de decisión.
P = length(NORMAL)+length(DTA);
labels= ones(P,1);  labels(NORMAL)= -1;

figure()
plot3(XS(NORMAL,1),XS(NORMAL,2),XS(NORMAL,3),'b*'); hold on;
h= plot3(XS(DTA,1),XS(DTA,2),XS(DTA,3),'rs'); hold off;
grid;
hAxis = get(h,'parent');
lims = axis(hAxis);

N=100;
[Xc, Yc, Zc]= meshgrid(linspace(lims(1),lims(2),N),linspace(lims(3),lims(4),N),linspace(lims(5),lims(6),N));

output= zeros(N,N,N);

tr_data= XS(:,[ 1 2 3 ]);
train= 1:P;
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','linear');
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','quadratic');
%svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function','polynomial');
svmStruct= svmtrain(tr_data(train,:), labels(train),'Kernel_Function',kernel_function);

for xc=1:N
    for yc=1:N
        for zc=1:N
            output(xc,yc,zc)= eval_svmStruct(svmStruct,[Xc(xc,yc,zc) Yc(xc,yc,zc) Zc(xc,yc,zc)]);
        end
    end
end

plot3(XS(NORMAL,1),XS(NORMAL,2),XS(NORMAL,3),'bo','MarkerSize',7); hold on;
h= plot3(XS(DTA,1),XS(DTA,2),XS(DTA,3),'rs','MarkerSize',7);
hAxis = get(h,'parent');
lims = axis(hAxis);

% Representamos los vectores de soporte.
hold on;
sv = svmStruct.SupportVectors;
% see if we need to unscale the data
scaleData = svmStruct.ScaleData;
if ~isempty(scaleData)
    for c = 1:size(sv, 2)
        sv(:,c) = (sv(:,c)./scaleData.scaleFactor(c)) - scaleData.shift(c);
    end
end
% plot support vectors
hSV = plot3(sv(:,1),sv(:,2),sv(:,3),'*k');

hpatch = patch(isosurface(Xc,Yc,Zc,output,0));
isonormals(Xc,Yc,Zc,output,hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
view([-126,36]);
axis tight
camlight left; 
set(gcf,'Renderer','zbuffer'); 
lighting phong;
grid;
hold off;

legend('NOR','AD','Support vectors','Decision surface');
title(['#SVs= ' num2str(size(sv,1))]);
y= 1;