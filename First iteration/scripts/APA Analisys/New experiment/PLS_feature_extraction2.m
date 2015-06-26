function XS= PLS_feature_extraction2(Y,X,ncomp)

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
