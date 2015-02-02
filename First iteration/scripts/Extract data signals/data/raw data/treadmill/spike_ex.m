function [mark]=spike_ex(mark,thres)
% [mark]=spike_ex(mark,thres)
% remove spikes and NaN from a marker vektor and interpolate
% data must be in columns
% does not work when last value is NaN
% can deal only with 2 dim arrays
% global m ma1 ma2 % ind2

ma1=isnan(mark);                              %  find NaN
% markold=mark;                                 % graphische Kontrolle, s.u.

s=size(mark);
ma2=zeros(s);
ma2(2:end,:)=abs(diff(mark))>thres;


m = ma1 | ma2;                                 % matrix where all outliers are 1 else 0

for col=1:s(2);
    if any(m(:,col));
        row=2;
        while row < s(1);
            if m(row,col)==1;                                     % first outlier detected
               anz=1;                                             % number of connected outliers
               while m(row+anz,col)==1 && row+anz+1 < s(1);        % check for end of array
                 anz=anz+1;
               end
%                if anz>1                                           % last difference is not an outlier
%                    anz=anz-1;
%                end
               a=(mark(row+anz,col)-mark(row-1,col))/(anz+1);     % increment
               new = (1:anz) * a + mark(row-1,col);               % interpolate new line
               anf=row+(col-1)*s(1);
               mark(anf:anf+anz-1)= new;
               row=row+anz;
            end
            row=row+1;
        end
    end
end

% figure(7)
% for col=1:s(2)
%     subplot(2,1,1)
%     plot(markold(:,col))
%     subplot(2,1,2)
%     plot(mark(:,col))
%     [d1,d2]=ginput(1);
% end
% %     