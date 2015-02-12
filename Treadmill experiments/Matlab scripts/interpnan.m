function data=interpnan(data);
% function data=interpnan(data);
% replace NaN by linear interpolation assuming equal sampling
%% doit
n=length(data);
idx=find(isnan(data));
while ~isempty(idx),
    idxgaps=find(diff(idx)>1,1);
    i1=idx(1)-1;
    if i1==0,
        d1=0;
        i1=1;
    else
        d1=data(i1);
    end
    if isempty(idxgaps),
        i2=idx(end)+1;
        if i2>length(data),
            d2=d1; % last data point, set to last known value
        else
            d2=data(i2);
        end
    else
        i2=idxgaps(1)+idx(1);
        d2=data(i2);
    end
    data(i1:i2)=d1+(d2-d1)*(0:i2-i1)/(i2-i1);
    idx=find(isnan(data));
end
data=data(1:n);
