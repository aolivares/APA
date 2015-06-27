function [ output ]= eval_svmStruct(svmStruct,sample)

if ~isempty(svmStruct.ScaleData)
    for c = 1:size(sample, 2)
        sample(:,c) = svmStruct.ScaleData.scaleFactor(c)*(sample(:,c) +  svmStruct.ScaleData.shift(c));
    end
end

sv = svmStruct.SupportVectors;
alphaHat = svmStruct.Alpha;
bias = svmStruct.Bias;
kfun = svmStruct.KernelFunction;
kfunargs = svmStruct.KernelFunctionArgs;

output = (feval(kfun,sv,sample,kfunargs{:})'*alphaHat(:)) + bias;
