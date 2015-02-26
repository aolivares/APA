% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||||| GAIT ANALYSIS ROUTINE ||||||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% read gaitWatch file plus data made by routine 'main' i.e. kalman-filtered
% pitch data of leg
% and then ... analyse gait!
% analyse soles data

clear
% 
% filenames ={'GaitWatch_1358_070214_1255.mat','insoles_subj12_01.csv';  % Martin Bartels
%             'GaitWatch_1359_070214_1257.mat','insoles_subj12_02.csv';
%             'GaitWatch_1360_070214_1300.mat','insoles_subj12_03.csv'};
        
        
% filenames ={'GaitWatch_1361_070214_1325.mat','insoles_subj13_01.csv';   % Soles data: Gyro& acc missing !! Subject Matthias Ertl
%             'GaitWatch_1362_070214_1327.mat','insoles_subj13_02.csv';   % Soles data : Gyro& acc missing
%             'GaitWatch_1363_070214_1329.mat','insoles_subj13_03.csv';   % soles data file missing
%             'GaitWatch_1364_070214_1333.mat','insoles_subj13_04.csv'};  % this works !!
%         
%         
% filenames ={'GaitWatch_1365_070214_1400.mat','insoles_subj_11_04.csv';   % Soles data : Gyro& acc missing Subject Lei Zang 
%             'GaitWatch_1366_070214_1402.mat','insoles_subj_11_05.csv';   % Soles data : Gyro& acc missing
%             'GaitWatch_1367_070214_1404.mat','insoles_subj_11_06.csv'};  % Soles data : Gyro& acc missing
% 
  filenames ={'GaitWatch_1320_040214_2016.mat','gw02.csv';   % This pair of files doesnt fit well !! together Pedro Sá
             'GaitWatch_1321_040214_2018.mat','gw03.csv';    % this seems to fit
             'GaitWatch_1322_040214_2022.mat','gw04.csv'};   % this pair seems to fit
        
         
         
vers = 2;    % take first data pair from this list



%  filename_soles =['C:\A-DATA\Qualisis&Treadmill\7th february treadmill\',filenames{vers,2}];
filename_soles =['C:\KB_DATA\Qualisis&Treadmill\4th february treadmill\ECnsoles dataset\',filenames{vers,2}];




filepath = 'C:\KB_DATA\GaitWatchData';
%[filename, filepath] = uigetfile('C:\A-DATA\GaitWatchData\*.mat', 'File auswählen');
filename = filenames{vers,1};
% load([filepath,'\',filename(1:end-4),'_KF.mat'])                                       % load Kalman filtered data ('KF_leg_pitch','acc_center','gy_center')
load(fullfile(filepath,filename));                                                     % load original data

% %  *********************   KF-Data Format   *******************************
% % kalman filter leg data (right shank right thigh, left shank, left thigh)
% KF_leg_pitch = [pitch_KF_right_shank, pitch_KF_right_thigh, pitch_KF_left_shank, pitch_KF_left_thigh];
% % acceleration trunk
% acc_center =[a_X_center_trunk_3_C', a_Y_center_trunk_3_C', a_Z_center_trunk_3_C'];
% % gyro trunk
% gy_center = [g_X_center_trunk_1_C,g_Y_center_trunk_1_C,g_Z_center_trunk_1_C];
% % gyro arm
% gy_arm = [g_X_right_arm_1_C, g_Y_right_arm_1_C, g_X_left_arm_1_C, g_Y_left_arm_1_C];

% 1. gY_r_shank x
% 2. aZ_r_shank X
% 3. aX_r_shank X
% 4. gY_r_thigh
% 5. aZ_r_thigh
% 6. aX_r_thigh
% 7. gY_l_shank x
% 8. aZ_l_shank X
% 9. aX_l_shank X
% 10.gY_l_thigh
% 11.aZ_l_thigh
% 12.aX_l_thigh
% 13.gY_l_arm  
% 14.gX_l_arm  
% 15.gY_r_arm  
% 16.gX_r_arm  
% 17.gZ_c_trunk
% 18.gY_c_trunk
% 19.gX_c_trunk
% 20.aZ_c_trunk
% 21.aY_c_trunk
% 22.aX_c_trunk
% 23.compass   


% ****************             calibrate GaitWatch data           % *******************

SampRate = double(FileHeader.SampFreq);
len_GW_data=length(data);

load CalibrationList.mat                                                          % calibrate data 
data = double(data);
lower_freq_limit=0.2;                                                             % highpass-Filter to compensate drift
[b2,a2]=butter(3,lower_freq_limit*2/(SampRate),'high');


for n=1:16
    if  FileHeader.ChannelNames(n,1) == 'a'                                       % accelerometers and gyros are calibrated differently according to alberto!
        data(:,n) = data(:,n) * k_array(n)  + bias_array(n);
    else
        data(:,n) = (data(:,n) - bias_array(n)) /  k_array(n)  ;
        data(:,n) = filtfilt(b2,a2,data(:,n));                                     % gyroscopes are high pass filtered to eliminate offset
    end
end

% ****************             keep only some GaitWatch data / create timeseries          % *******************

GW_data_index = [1:3 7:9];                                                         % GaitWatch data (only shank data, gyro and accelerometer)
GW_data_labels = FileHeader.ChannelNames(GW_data_index,:);
GW_data = timeseries(data(:,GW_data_index),(0:len_GW_data-1) / SampRate);

clear data


% ****************             read soles data         *******************
 
fid=fopen(filename_soles);
fmt='%*s %s%s%s%s %*s%*s%*s%*s%*s %s %*s%*s%*s %s%s%s%s %*s%*s%*s%*s%*s %s %*s%*s%*s';
soles_data_labels = textscan(fid, fmt,1,'Delimiter',',');                                 % header lesen 1 Zeile
fmt='%*d %d%d%d%d %*d%*d%*d%*d%*d %d %*d%*d%*d %d%d%d%d %*d%*d%*d%*d%*d %d %*d%*d%*d';
temp = textscan(fid,fmt,'Delimiter',',');                                           % header lesen 1 Zeile
fclose(fid);

soles_data = zeros(length(temp{1}),10);
for n=1:10
    soles_data(:,n) = temp{n};
end

clear temp

len_sdata=length(soles_data);
time_so=(0:len_sdata-1) / 20;                                                   % time in s / sampling rate = 20 Hz




% ****************             synchronize data by moving them manually          *******************

gd1 = GW_data.Data(:,[2 5]) - repmat(    mean(GW_data.Data(:,[2 5])),   length(GW_data.Data),   1);                        % Gyro-data aZ_r_shank aZ_l_shank 
gd1 = gd1 ./ repmat(max(gd1),length(gd1),1);

sd1 = soles_data(:,[5 10]) - repmat(     mean(soles_data(:,[5 10])), length(soles_data), 1);                                % Soles-data ADZ AIZ
sd1 = sd1 ./ repmat(max(sd1),length(sd1),1);


inc = 2000;                                                                   % display inc in Gyro data index
lag=0;
an=201;
en=an+inc-1;
taste = 1;




scrsz = get(groot,'ScreenSize');
figure(1)
set(gcf,'Position',[10 scrsz(4)/3 scrsz(3)-20 scrsz(4)/2]);


while an < len_GW_data
    hold off
    if en>len_GW_data
        en = len_GW_data;
    end
    [~,an_so] = min(abs((time_so)-GW_data.Time(an)));
    [~,en_so] = min(abs((time_so)-GW_data.Time(en)));
    if en_so > len_sdata, en_so = len_sdata; end
    clf
    for n=1:2
        subplot(2,1,n)
        plot(GW_data.Time(an:en),gd1(an:en,n),'linewidth',1.5)
        hold on
        plot(time_so(an_so:en_so),sd1(an_so:en_so,n),'linewidth',1.5)
        title('left button:left         right button: right        middle button:forward')
        axis([GW_data.Time(an) GW_data.Time(en) -1 1])
        legend({strrep(GW_data_labels(n+1,:),'_',' '), char(soles_data_labels{n*5})})
    end
       
    [~,~,taste] = ginput(1)
    switch taste
        case 3
            lag = lag + 0.05
            time_so = time_so + 0.05;
        case 1
            lag = lag - 0.05
            time_so = time_so - 0.05;
        case 2
            an=an+inc;
            en=an+inc-1;
    end
    
end

clear gd1 sd1

% ****************             Synchronize data           *******************


% now the time scalars of both data sets are synchronized
% construct 2 time series sets and synchronize data
% synchronize: 'union': time values are new and will be the same for both series
% Time of the gaitWatch-data will not be changed (but the indices of the GaitWatch data change because data are truncated)


 
soles_data = timeseries(soles_data,time_so);
[GW_data,soles_data] = synchronize(GW_data,soles_data,'union');   



% ****************             detect steps by looking at rotation of shank           *******************

clf

upper_freq_limit=4;                                                                                  % lowpass-Filter for peak detektion of shank-gyro data
[b,a]=butter(3,upper_freq_limit*2/(SampRate)); 

gd=[1 4];                                                                       % gyro-data
PkLocs=cell(2,1);                                                                                      % peak detektion shank gyros right and left leg (1,2)
for n=1:2
    subplot(2,1,n)
    plot(GW_data.Time, GW_data.Data(:,gd(n)))
    temp = filtfilt(b,a,GW_data.Data(:,gd(n)));
    [pks,locs] = findpeaks(temp,'minpeakheight',40);   
    weg=find(diff(locs) < SampRate/2)+1;                                                                % maximum 2 stps/s
    locs(weg)=[];
    PkLocs{n}=locs;                                                                                     % PkLocs = Indices der peaks
    text(GW_data.Time(PkLocs{n}),GW_data.Data(PkLocs{n},gd(n)),'+','HorizontalAlignment','center')
    legend(strrep(GW_data_labels(gd(n),:),'_',' '))
    xlabel('Time [ s ]')
end


title('Mark begin and end of epoch for evaluation')
[an_en,y]=ginput(2);
[~,an_en(1)] = min(abs((GW_data.Time)-an_en(1)));
[~,an_en(2)] = min(abs((GW_data.Time)-an_en(2)));
PkLocs{1}(PkLocs{1} < an_en(1) | PkLocs{1} > an_en(2)) = [];                   % delete steps that shall not be evaluated
PkLocs{2}(PkLocs{2} < an_en(1) | PkLocs{2} > an_en(2)) = [];




% % soles data zero-line shows drift and spikes!! This can be corrected
clf
for n=[1:4 6:9]
   plot(soles_data.Data(:,n))
   xlabel('mark zero-line')
   legend(soles_data_labels{n})
   [x,y] = ginput(1);
   ind = find(soles_data.Data(:,n) < y);
   soles_data.Data(ind,n) = y;
end

set(gcf,'Position',[200 50 scrsz(3)/2 scrsz(4)-110]);



% ****************             plot every step to see how it looks like           *******************
% mark IC and TC manually
% compute the descriptors of TC and IC in the filtered curve of the gyroscope
% PkLocs{1} right leg; PkLocs{2} left leg
%             index 2: 5 instances on the filtered gyro-curve:
%                      peak/minimum before peak/TC/minimum after peak/IC
%             index 3: index / value

PkLocs{1}(:,5,2) = 0;
PkLocs{2}(:,5,2) = 0;

an_en(1) = 0.5*SampRate;        % No of values before peak
an_en(2) = 0.5*SampRate;        % No of values after peak


for side=1:2
    temp = filtfilt(b,a,GW_data.Data(:,gd(side)));
    for step=PkLocs{side}(:,1,1)'
        an = step - an_en(1);
        en = step + an_en(2);
        subplot(3,1,1)                                       % plot acceleration
        if side == 1
            plot(an:en, GW_data.Data(an:en,[2 3]),'linewidth',1.5)
            legend(strrep(GW_data_labels(2,:),'_',' '),strrep(GW_data_labels(3,:),'_',' '))
        else
            plot(an:en, GW_data.Data(an:en,[5 6]),'linewidth',1.5)
            legend(strrep(GW_data_labels(5,:),'_',' '),strrep(GW_data_labels(6,:),'_',' '))
        end
        
        
        subplot(3,1,2)                                       % plot Gyro
            hold off
            plot(an:en,GW_data.Data(an:en,gd(side)),'linewidth',1.5)
            hold on
            plot(an:en,temp(an:en),'linewidth',1.5)
            legend(strrep(GW_data_labels(gd(side),:),'_',' '),'LP-Filt 4 Hz')
               
        subplot(3,1,3)                                       % plot soles
        if side == 1
            plot(an:en, soles_data.Data(an:en,1:4),'linewidth',1.5)
            legend(char(soles_data_labels{1}), char(soles_data_labels{2}), char(soles_data_labels{3}), char(soles_data_labels{4}))
        else
            plot(an:en, soles_data.Data(an:en,6:9),'linewidth',1.5)
            legend(char(soles_data_labels{6}), char(soles_data_labels{7}), char(soles_data_labels{8}), char(soles_data_labels{9}))
        end
        xlabel('mark TC and IC')
       
        % find the descriptors of the filtered gyro curve
        ind = find(PkLocs{side}(:,1,1) == step);                              % find index 
        PkLocs{side}(ind,1,2) = temp(step);                                     % value of peak
       
        [value,index] = min(temp(an:step));                                    % find Minimum before peak
        PkLocs{side}(ind,2,1) = an + index - 1;                                 % index of minimum bfore peak
        PkLocs{side}(ind,2,2) = value;                                          % value of minimum before peak
        
        [value,index] = min(temp(step:en));                                    % find Minimum after peak
        PkLocs{side}(ind,4,1) = step + index - 1;                               % index of minimum after peak
        PkLocs{side}(ind,4,2) = value;                                          % value of minimum after peak
               
       [x,~] = ginput(2);
            
        PkLocs{side}(ind,3,1) = round(x(1));                               % index of terminal contact
        PkLocs{side}(ind,3,2) = temp(round(x(1)));;                        % value of terminal contact
             
        PkLocs{side}(ind,5,1) = round(x(2));                               % index of initial contact
        PkLocs{side}(ind,5,2) = temp(round(x(2)));;                        % value of initial contact
   
     end
end


% C:\KB_Data\Qualisis&Treadmill\Evaluation\Soles\
% save(['C:\KB_Data\Qualisis&Treadmill\Evaluation\Soles\',filename(1:end-4),'_GW_Soles.mat'],...
%             'filename',...
%             'filename_soles',...
%             'GW_data','soles_data','PkLocs')
        

% figure
% clf
% subplot(2,1,1)
% plot(GW_data.Data(1:2000,[2 3]),'linewidth',1.5)
% legend(strrep(GW_data_labels(2,:),'_',' '),strrep(GW_data_labels(3,:),'_',' '))
% subplot(2,1,2)
% plot(GW_data.Data(1:2000,[5 6]),'linewidth',1.5)
% legend(strrep(GW_data_labels(5,:),'_',' '),strrep(GW_data_labels(6,:),'_',' '))



