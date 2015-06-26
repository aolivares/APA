function [acc, sen, spe] = SVM_LOU(data,labels)

P = size(data,1);
N = P;  % Número de veces que se realiza el experimento Leave-M-Out = numero de grupos disjuntos

NPos = sum(labels==1);
NNeg = sum(labels~=1);


kern = {'linear','quadratic','polynomial','rbf'};

for k = 1:length(kern)
    kernel = char(kern(k));
    Errors = ones(N,1);
    for n = 1:N
        
        
        train_set = true(1,P); train_set(n) = false; test_set =~train_set;
        train_labels = labels(train_set); test_labels = labels(test_set);
        
        tr_data = data(train_set,:);
        te_data = data(test_set,:);
        %         [train_data, test_data] = sort_pca_FDR (train_data, test_data, tr_labels(train_set),PC);
        
        
        %             svmStruct = svmtrain(tr_data, train_labels,'Kernel_Function', tp,'RBF_Sigma', Sigma);
        svmStruct = svmtrain(tr_data, train_labels, 'Kernel_Function', kernel,'Showplot', false);
        class = svmclassify(svmStruct, te_data)';
        
        
        % Leave-M-Out
        class = cast(class,'double');
        Errors(n) = length(find(class ~= test_labels'));
        
        %         if Errors(n,f) == 1
        %             if te_labels
        %                 FN(f) = FN+1;
        %             else
        %                 FP = PF+1;
        %             end
        %         else
        %             if te_labels
        %                 TP = TP+1;
        %             else
        %                 TN = TN+1;
        %             end
        %         end
        
    end
    
    
    % Leave-M-Out
    Total_Errors = sum(Errors,1);
    Total_test = N;
    acc(k) = (Total_test - Total_Errors)./Total_test * 100;
    
    POS = labels==1;
    NEG = labels~=1;
    FP = sum(Errors(POS,:),1); TP = NPos-FP;
    FN = sum(Errors(NEG,:),1); TN = NNeg-FN;
    sen(k) = (TP./(TP+FN));
    spe(k) = (TN./(TN+FP));
    
end





