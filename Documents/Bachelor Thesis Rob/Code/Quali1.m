% Quali1.m
% Plot Qualisis data


clear

% fn = 'gw_subj11_05';                                                  % filename
fn = 'gw_subj11_06';                                                  % filename
fp = 'C:\\A-DATA\\Qualisis&Treadmill\\7th february treadmill';        % filepathe

load (sprintf('%s\\%s.mat',fp,fn))

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

x=1;z=3;m1=1;m2=2;m3=3;

% right shank
seg_1 = atan2(squeeze(data(m1,x,:)-data(m2,x,:)), squeeze(data(m1,z,:)-data(m2,z,:)));        % marker 1 vs Marker 2 x/z


% seg_1 = atan(squeeze((data(m1,x,:)-data(m2,x,:))./(data(m1,z,:)-data(m2,z,:))));        % marker 1 vs Marker 2 x/z
seg_2 = atan(squeeze(       (data(m3,z,:)-data(m2,z,:)) ./ (data(m2,x,:)-data(m3,x,:))  ));        % marker 2 vs Marker 3 x/z
seg_3 = atan(squeeze(       (data(m1,z,:)-data(m3,z,:)) ./ (data(m3,x,:)-data(m1,x,:))  ));        % marker 3 vs Marker 1 


Q_leg_pitch = mean(seg_1 + seg_2 + seg_3);

seg_1b = atan2(squeeze(data(m1,x,:)-data(m2,x,:)), squeeze(data(m1,z,:)-data(m2,z,:)));        % marker 1 vs Marker 2 x/z

figure
% plot(Q_leg_pitch)
plot(seg_1)



figure
plot(seg_1b)

hold on
plot(seg_2,'r')
plot(seg_3,'k')





















% construct axes

clf
mi = min(min(shiftdim(data,2)));   % maximum x-y-z
ma = max(max(shiftdim(data,2)));   % minimum x-y-z
axis([mi(1) ma(1) mi(2) ma(2) mi(3) ma(3)])
set(gca,'ZDir','reverse','xgrid','on','ygrid','on','zgrid','on')
axis equal

xlabel('X')
ylabel('Y')
zlabel('Z')
hold on

% show qualisys markers in their initial position

size_marker = 20;               % radius of graphic symbol for marker
[X,Y,Z] = sphere(30);           % graphic symbol for plotting the markers
X = X * size_marker;
Y = Y * size_marker;
Z = Z * size_marker;
% b=surf(X,Y,Z,c,'LineStyle','none','FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5);


c = ones(size(X));
colormap([1  0  0])             % red markers

marker=cell(n_marker,1);                % graphic object marker (sphere)
for n=1:n_marker                        % bring 15 markers into their initial position
    x = X + data(n,1,1);                
    y = Y + data(n,2,1);
    z = Z + data(n,3,1);
    marker{n} = surf(x,y,z,c,'LineStyle','none','FaceLighting','phong','FaceColor',[1 0 0],'AmbientStrength',0.5);
end
light('Position',[1000 -200 -5000],'Style','infinite');




% show 3D movie (this is very very slow!!)


for n2 = 1500:10:n_frames
    for n=1:n_marker
        x = X + data(n,1,n2);
        y = Y + data(n,2,n2);
        z = Z + data(n,3,n2);
        set(marker{n},'XData',x)
        set(marker{n},'YData',y)
        set(marker{n},'ZData',z)
        drawnow
    end
end


% for to check how well the markers were fixed to the body, now for each marker-triangel, the distances of the markers is calculated 
% It is assumed that this will not vary much if they were well fixed.
% this would be a quality control

dist_var = zeros(n_marker,1);                                                  % variance of the distance-data


figure
farbe ={'r','g','b'};
clf
ind = 0;                                     %  index
ind2 = 0;
% ind3 = 1;                                    % index for variance-array
for n=1:5                                    % 5 triangles
    subplot(5,1,n)
    for n2 =1:3                             % each has 3 distances
        ind = ind + 1;
%         ind3 = ind3 + 1;
        if n2 < 3, ind2 = ind + 1; else ind2 = ind - 2; end
        dist = squeeze(data(ind,:,:)) - squeeze(data(ind2,:,:));        % subtract position of markers to get the difference vector
        dist = sqrt(dot(dist,dist));                                    % compute the norm (length) of each difference vector
        plot(dist-mean(dist),farbe{n2})
        dist_var(ind) = var(dist);
        hold on
    end
    title([char(data_labels{ind-2}),'     ',char(data_labels{ind-1}),'     ',char(data_labels{ind})])
%     ind3 = ind3 + 1;
end

axes('position',[ 0 0 1 1], 'visible','off')
text(0.5,0.97,'Distance between markers of one triangle','units','normalized','FontSize',12,'FontWeight','bold','HorizontalAlignment','center')



% show variance of the data in a barplot
figure
subplot(2,1,1)
bar(dist_var)
title('Variance of distance of Qualisis-markers')
ylabel('Variance in mm')

axes('position',[ 0 0 1 1], 'visible','off')             % show the name of the markers

for ind=1:8
    if mod(ind,3) ~= 0, ind2 = ind + 1; else ind2 = ind - 2; end
    disp([ind ind2])
    text(0.125,0.5 - ind * 0.05,sprintf('%d.) %s - %s', ind, data_labels{ind},data_labels{ind2}),...
        'units','normalized')
 end
for ind=9:15
    if mod(ind,3) ~= 0, ind2 = ind + 1; else ind2 = ind - 2; end
    disp([ind ind2])
    text(0.6,0.475 - (ind-8) * 0.05,sprintf('%d.) %s - %s', ind, data_labels{ind},data_labels{ind2}),...
        'units','normalized')
end    
  

 
% for to check how well the markers were fixed to the body, compute the frequency spectrum of each marker
% It is assumed that high frequencies will be cused by marker not being loose 

figure
farbe ={'r','g','b'};
clf
ind = 0;                                     %  index
lower_freq_limit=0.5;  
[b,a]=butter(4,lower_freq_limit*2/(SampRate),'high');

        
for n=1:5                                      % 5 triangles
    subplot(3,2,n)
    for n2 =1:3                                % each triangle has 3 markers
        ind = ind + 1;
        dist = squeeze(data(ind,:,:));                              % x y z of one marker
        dist = dist - repmat(mean(dist')',1,n_frames);              % get rid of the mean
        dist = dot(dist,dist);                                      % square and sum up
        dist = filtfilt(b,a,dist);                                  % highpass filter
        [Pxx,F] = pwelch(dist,[],[],[],SampRate);                   % [Pxx,F] = pwelch(X,WINDOW,NOVERLAP,NFFT,Fs)
        plot(F(1:100),Pxx(1:100),farbe{n2})
        hold on
    end
    % axis([0 F(100) 0 10.^9]) 
    legend(char(data_labels{ind-2}),char(data_labels{ind-1}),char(data_labels{ind}))
end
axes('position',[ 0 0 1 1], 'visible','off')
text(0.5,0.97,'Frequency spectra of individual markers','units','normalized','FontSize',12,'FontWeight','bold','HorizontalAlignment','center')
       






% calculate position of hip joint and knee joint first in 2d space
% set initial position of marker 5 as a origin and orientation of marker 5 - 4 as veritcal
% plot all markers of right leg and hip in the coordinate system of 4/5/6 
% hip and shank will rotate in this coordinate system
% determine rotation axis by aproximation
% do this for 101 frames


% step 1: select epoch for calculation##
figure

for mk=7:9
    subplot(3,1,mk-6)
   plot(squeeze(data(mk,1:3,:))')
    title(data_labels{mk})
end
[start_ind,dummy]=ginput(1);
start_ind = fix(start_ind);



clf
mi = min(min(shiftdim(data,2)));   % maximum x-y-z
ma = max(max(shiftdim(data,2)));   % minimum x-y-z
% axis([mi(1) ma(1) mi(2) ma(2) mi(3) ma(3)])
axis([-1000 0 -600 0 -1200 0])
set(gca,'ZDir','reverse','xgrid','on','ygrid','on','zgrid','on')
axis equal

xlabel('X')
ylabel('Y')
zlabel('Z')
hold on

% show qualisys markers in their initial position

size_marker = 20;               % radius of graphic symbol for marker
[X,Y,Z] = sphere(30);           % graphic symbol for plotting the markers
X = X * size_marker;
Y = Y * size_marker;
Z = Z * size_marker;
c = ones(size(X));
colormap([1  0  0])             % red markers
% b=surf(X,Y,Z,c,'LineStyle','none','FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5);


marker=cell(n_marker,1);                % graphic object marker (sphere)
for n=1:9                        % bring 15 markers into their initial position
    x = X + data(n,1,1);                
    y = Y + data(n,2,1);
    z = Z + data(n,3,1);
    marker{n} = surf(x,y,z,c,'LineStyle','none','FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5);
end
light('Position',[1000 -200 -5000],'Style','infinite');




new_data = zeros(n_marker,4,400);                                                               % rotated data
% start_ind = 4000;

for ind = 1:400                                                                                 % ind belongs to new_data
    n = ind*10+start_ind;                                                                       % n belongs to data
    TM1 = makehgtform('translate',-data(5,1,n),-data(5,2,n),-data(5,3,n));                      % first step: calculate translation matrix to set marker 5 into center
    alpha = asin(  (data(5,1,n) - data(4,1,n)) /  norm( data(5,[1 3],n)'-data(4,[1 3],n)' ) );  % second step: determine rotation around y  (sin alpha = dx/hy)
    beta = -asin(  (data(5,2,n) - data(4,2,n)) /  norm( data(5,[2 3],n)'-data(4,[2 3],n)' ) );   % third step: determine rotation around x  (sin beta = dy/hy)
    gamma = -asin(  (data(5,2,n) - data(6,2,n)) /  norm( data(5,[1 2],n)'-data(6,[1 2],n)' ) );   % fourth step: determine rotation around z  (sin gamma = dy/hy)
    fprintf('alpha : %6.2f Grad;   beta : %6.2f Grad;   gamma : %6.2f Grad\n',alpha / pi * 180, beta / pi * 180, gamma / pi * 180)
    RM1 = makehgtform('yrotate',alpha,'xrotate',beta,'zrotate',gamma);    % ,'zrotate',gamma 
    TM2 = makehgtform('translate',data(5,1,start_ind),data(5,2,start_ind),data(5,3,start_ind));   % last step: calculate translation matrix to set marker 5 into position of initial marker 5
    for mk=1:9
        new_data(mk,:,ind) = TM2 * RM1 * TM1 * data(mk,:,n)';                  % compute new data
        set(marker{mk},'XData',X + new_data(mk,1,ind),'YData',Y + new_data(mk,2,ind),'ZData',Z + new_data(mk,3,ind))
    end
    drawnow
end
    
 
% compute the center of rotation of the hip markers
% we have 3 hip markers and we do this only in the x-z plane
clf
right_hip_coords=zeros(3,3);      % (mk,x/z/R)

col = {'r','g','b'};  
figure

for mk=1:3                          % marker 7-9!!
   [right_hip_coords(mk,1),right_hip_coords(mk,2),right_hip_coords(mk,3),a] = circfit(new_data(mk+6,1,:),new_data(mk+6,3,:));
   subplot(1,3,1)
   plot(squeeze(new_data(mk+6,1,:)),squeeze(new_data(mk+6,3,:)),['.',col{mk}])
   hold on
   plot(squeeze(new_data(mk+3,1,:)),squeeze(new_data(mk+3,3,:)),['.',col{mk}])
   viscircles([right_hip_coords(mk,1) right_hip_coords(mk,2)],right_hip_coords(mk,3),'EdgeColor',col{mk})
   
   
   subplot(1,3,2)
   plot(squeeze(new_data(mk+6,2,:)),squeeze(new_data(mk+6,3,:)),['.',col{mk}])
   hold on
   plot(squeeze(new_data(mk+3,2,:)),squeeze(new_data(mk+3,3,:)),['.',col{mk}])
   
   subplot(1,3,3)
   plot(squeeze(new_data(mk+6,2,:)),squeeze(new_data(mk+6,1,:)),['.',col{mk}])
   hold on
   
end
 subplot(1,3,1)
 xlabel('X')
 ylabel('Z')
 set(gca,'YDir','reverse'), axis equal                 
                
subplot(1,3,2)
 xlabel('Y')
 ylabel('Z')
 set(gca,'YDir','reverse'), axis equal                 
                                
 subplot(1,3,3)
 xlabel('Y')
 ylabel('X')
 axis equal    
                
              

[xc,yc,R,a] = circfit(new_data(9,1,:),new_data(9,3,:));
d = ones(size(X))+1;
hip_joint = surf(X+xc,Y-450,Z+yc,'LineStyle','none','FaceLighting','phong','FaceColor',[0 1 0],'AmbientStrength',0.5);










