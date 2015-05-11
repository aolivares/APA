% -------------------------------------------------------------------------
% 1) Load original GaitWatch data
% and Kalman filtered GaitWatch data which were produced by progarm main.m
% -------------------------------------------------------------------------

load('GaitWatch_data.mat')

KF_leg_pitch = [pitch_KF_right_shank, pitch_KF_right_thigh, pitch_KF_left_shank, pitch_KF_left_thigh];

% -------------------------------------------------------------------------
% 1) Load corresponding Qualisys data
% -------------------------------------------------------------------------

fn = 'gw_subj5_04';
load('gw_subj5_04.mat')

% eval([fn,'.Trajectories.Labeled.Labels{1}']);
SampRate = eval([fn,'.FrameRate']);
data = eval([fn,'.Trajectories.Labeled.Data']);
% data(:,4,:) = [];                                                      % skip trace no 4, it doesnt contain information
data = data * -1;                                                      % x, y and z must be inverted to bring it into our coordinate system
data(:,4,:) = 1;                                                         % set all values of dim 4 to 1 to allow for homogeneous matrix multiplications 
[n_marker,dummy,n_frames] = size(data);
data_labels = eval([fn,'.Trajectories.Labeled.Labels']);
for n=1:n_marker, data_labels{n} = strrep(data_labels{n},'_',' '); end                  % remove _ because this makes trouble when printed
eval(['clear ' fn]);                                                   % remove struct for memory reasons




% show how many data are missing in each channel

fprintf('\n\n')
for n=1:n_marker
    fprintf('Marker %d.) %s : No of NaN : %d \n', n, data_labels{n}, sum(isnan(squeeze(data(n,1,:)))))
end
fprintf('\n\n')


% plot raw data and interpolate missing values

mk = 1;

while mk
    mk = input('Do you want to see and correct raw data? Type marker number for yes or 0 for no : ');
    if mk
        clf
        for n=1:3
            subplot(3,1,n)
            if n ==  1, title(data_labels{mk}), end
            hold on
            d2 = interpnan(squeeze(data(mk,n,:)));                            % interpolated data in red
            plot(d2,'r')
            plot(squeeze(data(mk,n,:)),'b')                                   % old data in blue
        end
        ind = input('Shall the missing data be interpolated? (0/1)  : ');     % correct old data
        if ind 
            for n2 =1:3
                data(mk,n2,:) = interpnan(squeeze(data(mk,n2,:)));
            end
        end
    end
end
            




% calculate angular position of each marker trialgle in 2D
% and put them into data array Q_leg_pitch
% [pitch_KF_right_shank, pitch_KF_right_thigh, pitch_KF_left_shank, pitch_KF_left_thigh];

Q_leg_pitch = zeros(length(data),4);


order = {   'right lo shank','right up shank','right back shank';...
            'right lo thigh','right up thigh','right back thigh';...
            'left lo shank','left up shank','left back shank';...
            'left lo thigh','left up thigh','left back thigh';...
            'right hip','back hip','left hip'};
x=1;z=3;lag = zeros(4,1);
figure(1)
%  big_ax=axes('position', [0 0 1 1],'visible', 'off');
    
for n=1:4                                                                       % 4 leg segments
    m1=0;m2=0;m3=0;
    
    for n2=1:3                                                                  % 3 markers
        ind = find(strcmp(order{n,n2},data_labels));
        if isempty(ind)
            fprintf('Error: cannot find match for %s\n',order{n,n2})
            beep
        else
            switch n2
                case 1, m1 = ind;
                case 2, m2 = ind;
                case 3, m3 = ind;
            end
        end
    end
    
    
    
    seg_1 = -atan(squeeze((data(m2,x,:)-data(m1,x,:))./(data(m1,z,:)-data(m2,z,:))));        % marker 1,2: x/z
    seg_2 = -atan(squeeze((data(m2,z,:)-data(m3,z,:))./(data(m2,x,:)-data(m3,x,:))));        % marker 2,3: z/x
    seg_3 = -atan(squeeze((data(m1,z,:)-data(m3,z,:))./(data(m1,x,:)-data(m3,x,:))));        % marker 3,1: x/z

     Q_leg_pitch(:,n) = mean([seg_1 seg_2 seg_3],2); 
%     
%     
%     figure(1), clf
%     plot(  data([m1 m2 m3],1,1),   data([m1 m2 m3],3,1),'o')
%     text(data([m1 m2 m3],1,1),   data([m1 m2 m3],3,1), {'1','2','3'})
%     set(gca,'YDir','reverse','xgrid','on','ygrid','on')
%     axis equal
%     
    
    
    figure(1), clf
    subplot(4,1,1)
    
    plot(seg_1)
    title([data_labels{m1} '    ' data_labels{m2} '    ' data_labels{m3}]) 
    
    subplot(5,1,2)
    plot(seg_2)
    
    subplot(5,1,3)
    plot(seg_3)
    
    subplot(5,1,4)
    plot(Q_leg_pitch(:,n))
    title('Mean of 3 segments of one marker triangle')
    
    subplot(5,1,5)
    plot(KF_leg_pitch(:,n)) %KF_leg_pitch
    title('KF leg pitch')
    
    
    if n == 1
        axes('position', [0 0 1 1],'visible', 'off')
        t=text(0.5, 0.05,'Mark epoch of quiet standing in window 4: 2 clicks','HorizontalAlignment','center','FontSize',12);
        [an_en_q,~] = ginput(2);
        an_en_q = round(an_en_q);
        
        set(t, 'string','Mark epoch of quiet standing in window 5: 2 clicks')
        [an_en_kf,~] = ginput(2);
        an_en_kf = round(an_en_kf);
        
    else
        axes('position', [0 0 1 1],'visible', 'off')
        text(0.5, 0.05,'Click one time anywhere','HorizontalAlignment','center','FontSize',12);
        [~,~] = ginput(1);
    end
    
     Q_leg_pitch(:,n) = Q_leg_pitch(:,n) - mean(Q_leg_pitch(an_en_q(1):an_en_q(2),n));                            % position is set to zero during quiet standing
    KF_leg_pitch(:,n) = KF_leg_pitch(:,n) - mean(KF_leg_pitch(an_en_kf(1):an_en_kf(2),n));                       % position is set to zero during quiet standing
    
    
    
     maxlag=round(length(data)/5);
    [c,lags] = xcorr(Q_leg_pitch(:,n),KF_leg_pitch(:,n),maxlag);
    [~,ind] = max(c);
    lag(n) = lags(ind);
    
    
    figure(2),clf
    plot(KF_leg_pitch(-lag+1:end,n),'r')
    hold on
    plot(Q_leg_pitch(:,n)*(180/pi))
        
        
    [~,~] = ginput(1);    

    
    
end
        

lag=round(mean(lag));
if lag < 0
    KF_leg_pitch = KF_leg_pitch(-lag+1:end,:);                                                             % align both data
else
    Q_leg_pitch = Q_leg_pitch(-lag+1:end,:);                                    % falsch, index wird negativ, muss +lag heissen
end

if length(KF_leg_pitch) > length(Q_leg_pitch)                                                              % cut the longer of both
    KF_leg_pitch = KF_leg_pitch(1:length(Q_leg_pitch),:);
end

if length(Q_leg_pitch) > length(KF_leg_pitch)                                                              % both should have the same size now
    Q_leg_pitch = Q_leg_pitch(1:length(KF_leg_pitch),:);
end
    

% plot(KF_leg_pitch(-lag+1:end,n),'r')

% To do:
% arrays k?rzen und auf eine L?nge bringen
% Differenz berechnen
% Kennwerte der Marker notieren

for n=1:4
    clf
    subplot(2,1,1)
    plot(KF_leg_pitch(:,n),'r')
    hold on
    plot(Q_leg_pitch(:,n)*(180/pi))
    
    subplot(2,1,2)
    plot(KF_leg_pitch(:,n)-Q_leg_pitch(:,n)*(180/pi))
    [~,~] = ginput(1); 
end




pitch_QS_right_shank = Q_leg_pitch(:, 1)';
pitch_QS_right_thigh = Q_leg_pitch(:, 2)';
pitch_QS_left_shank = Q_leg_pitch(:, 3)';
pitch_QS_left_thigh = Q_leg_pitch(:, 4)';

save('Data/Qualisys_data', 'pitch_QS_right_shank', ...
     'pitch_QS_right_thigh', 'pitch_QS_left_shank', ...
     'pitch_QS_left_thigh');

