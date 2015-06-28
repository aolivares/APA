function [X,Y,auc]=roc(data,labels,step)

% Plot the ROC curve of the data using SVM classification
% Inputs:
%  - data  (double)  data input containing the feature vectors for SVM
%                    classficatons  with dimensions #observations times #features.
%  - labels (logical) vector containing the labels of the feature vectors 
%                     with the same dimensions as #observations, and 1 for Positive and 0 for negative
%  - step (double) number of total points in the ROC curve
%
% Outputs:
%  - X (double) vector containing the x axes values of the ROC
%               (1-specificity)
%  - Y (double) vector containing the y axes values of the ROC
%                (sensitivity)
%  - auc (double) Area Under the Curve (AUC) ROC, calculated using a
%                 rectangular aproximation

% Ignacio Alvarez 5/4/2011
%nfolds=103;




P=size(data,1);
kernels= {'linear' }; % 'rbf''quadratic'
output=zeros(P,1);
% indices = crossvalind('Kfold',P,nfolds);
for nk=1:numel(kernels)

    % Calculate the scores using the SVM
    for p=1:P
    test = false(P,1); test(p)=true;  train = ~test;
    %   test = (indices == p); train = ~test;
    sample=data(test,:);
    tt_labels=labels(train);
    svmStruct = svmtrain(data(train,:),labels(train),'Kernel_Function',char(kernels(nk)));
    if ~isempty(svmStruct.ScaleData)
    for c = 1:size(sample, 2)
        sample(:,c) = svmStruct.ScaleData.scaleFactor(c)*(sample(:,c) +  svmStruct.ScaleData.shift(c));
    end
    end
    output(p)=svmStruct.KernelFunction(svmStruct.SupportVectors,sample)'*svmStruct.Alpha+svmStruct.Bias;
    end

    
    % Define intial variables for step, counter and AUC
    max_val=max(output);
    min_val=min(output);
    cont=1;
    X=zeros(1,step+1);
    Y=zeros(1,step+1);
    % Calculate output classes compared to a varying score and obtain spe
    % and sen using classperf
    for i=min_val:(max_val-min_val)/step:max_val
        if all(svmStruct.Alpha>0==labels(svmStruct.SupportVectorIndices)')
        classes=int8(output>i)';
        elseif all(svmStruct.Alpha<0==labels(svmStruct.SupportVectorIndices)')
            classes=int8(output<i)';
        else
            disp('Raro!')
        end
        cp=classperf(labels,classes,'Positive',1,'Negative',0);
        X(cont)=1-cp.Specificity;
        Y(cont)=cp.Sensitivity;
        
        
        cont=cont+1;

    end
%     Xpar(p,:)=X;
%     Ypar(p,:)=Y;
end
%     X=mean(Xpar);
%     Y=mean(Ypar);
    % Calcula la contribución al área del rectaungulo de base step y de
        % altura sen
    
        auc=(X(2:end)-X(1:end-1))*Y(1:end-1)'+(X(2:end)-X(1:end-1))*(Y(2:end)-Y(1:end-1))'/2
    plot(X,Y)
    hold on
    
end
