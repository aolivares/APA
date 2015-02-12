% Quali1.m
% Plot Qualisis data


clear

fn = 'gw_subj11_05';                                                  % filename
fp = 'C:\Users\Usuario\Desktop\Experiments in Munich\7th february treadmill';        % filepathe

load (sprintf('%s\\%s.mat',fp,fn))

% eval([fn,'.Trajectories.Labeled.Labels{1}']);
SampRate = eval([fn,'.FrameRate']);
data = eval([fn,'.Trajectories.Labeled.Data']);
data(:,4,:) = [];                                                      % skip trace no 4, it doesnt contain information
data = data * -1;                                                      % x, y and z must be inverted to bring it into our coordinate system
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
            




% construct axes

clf
mi = min(min(shiftdim(data,2)));   % maximum x-y-z
ma = max(max(shiftdim(data,2)));   % minimum x-y-z
axis([mi(1) ma(1) mi(2) ma(2) mi(3) ma(3)])
set(gca,'ZDir','reverse')
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
b=surf(X,Y,Z,c,'LineStyle','none','FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5);


marker=cell(n_marker,1);                % graphic object marker (sphere)
for n=1:n_marker                        % bring 15 markers into their initial position
    x = X + data(n,1,1);                
    y = Y + data(n,2,1);
    z = Z + data(n,3,1);
    marker{n} = surf(x,y,z,c,'LineStyle','none','FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5);
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
       



