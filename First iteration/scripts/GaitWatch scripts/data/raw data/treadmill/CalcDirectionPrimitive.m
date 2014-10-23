function [mag_data,mag_time_ind] = CalcDirectionPrimitive(compass_trace)
% CalcDirectionPrimitive
% direction is in arcus not degree
% mag_time_ind = indexes of mag data, used for plotting them
% calculate direction of patient with GaitWatch.
% Does not take ino accuont pitch
% and is therefore not very precise

mag_cal =[340.6229  316.3078;                 % preliminary cal data and offset for x and y
         -425.7786 -226.0513;
                0         0];
mag=cell(3,1);                                 % only indices

offset1 = 20000;
offset2 = 32764;

mag{1} = find(compass_trace < -offset1/2);                                                   % Indices
mag{2} = find(compass_trace >= -offset1/2 & compass_trace < offset1/2);
mag{3} = find(compass_trace >= offset1/2 & compass_trace < offset2);

min_val = min([length(mag{1}), length(mag{2}), length(mag{3})]);
mag_time_ind = mag{1}(1:min_val);                                                   % 

mag_data = double([compass_trace(mag{1}(1:min_val))+ offset1,...
                   compass_trace(mag{2}(1:min_val)),...
                   compass_trace(mag{3}(1:min_val))- offset1]);
               
                             
for n=1:2                                                                                     % remove spikes by interpolation
    mag_data(:,n) = spike_ex_2(mag_data(:,n),100);
    mag_data(:,n) = (mag_data(:,n) - mag_cal(n,2)) / mag_cal(n,1);
end

mag_data(:,1) = atan2(mag_data(:,1),mag_data(:,2));
mag_data(2:end,2)=diff(mag_data(:,1));
mag_data(1,2)=0;
ind = find((mag_data(:,2) < -1) | (mag_data(:,2) > 1));                                    % kann aber auch grösser sein
mag_data(ind,2) = mag_data(ind-1,2);
mag_data=cumsum(mag_data(:,2));
end
