% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||| QUALISYS DATA PROCESSING ROUTINE |||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% The following file contains the routine which reads the data gathered
% using the Qualisys optical motion tracking system. It is prepared to
% process two different kinds of measurements: those having only two
% markers per segment, and those having three markers per segment. For the
% first case, it is only able to compute the pitch angle, while for the
% second case it computes both the pitch and roll angles. 

% -------------------------------------------------------------------------
% Authors: Alberto Olivares & Kai Bötzel.
% Entity: Universidad de Granada & Ludwig-Maximilians Universität München.
% Last modification: 13/12/2013.
% -------------------------------------------------------------------------

% 1) Definition of the flag which controls if the figures showing the X, Y
%    and Z coordinates of all the markers are shown.
showXYZ = 'yes';

% 2) The user is prompted to select the file he wishes to load.
[FileName,PathName] = uigetfile('*.mat','Select Qualisys data file');
structName = open(fullfile(PathName,FileName));
names = fieldnames(structName);

% 3) Extract the marker labels, the sampling frequency and the number of
%    total frames. 
eval(strcat('marker_labels=structName.',names{1},...
    '.Trajectories.Labeled.Labels;'));
eval(strcat('f=structName.',names{1},'.FrameRate;'));
eval(strcat('num_frames=structName.',names{1},'.Frames;'));

% 4) Build the time signal.
time_qs=zeros(num_frames,1);
for i=1:num_frames-1
 time_qs(i+1)=time_qs(i)+1/f;
end

% 5) Plot the X, Y and Z coordinates of each one of the markers.
if strcmpi(showXYZ,'yes')
    for i=1:length(marker_labels)
        figure
        plot(time_qs,squeeze(eval(strcat('structName.',names{1},...
            '.Trajectories.Labeled.Data(',num2str(i),',1:3,:)')))')
        title(marker_labels{i})
        legend('X axis','Y axis','Z axis')
        xlabel('Time (s)')
        ylabel('Position (mm)')
    end
end

% 6) Build a list showing all the markers. The user has to select which
%    markers he wants to process. There are two possible selections: 2 
%    markers or 3 markers.

ok_flag = 0;
while ok_flag == 0
    segment_select = listdlg('ListString',marker_labels,'Name',...
        'Select the markers','ListSize',[240 150],'SelectionMode',...
        'multiple');
    if length(segment_select) == 2 || length(segment_select) == 3 || ...
            length(segment_select) == 4
        ok_flag = 1;
    else
        msgbox(sprintf(['You selected %d markers. You must select ',...
            'either 2, 3 or 4 markers.'],length(segment_select)));
    end
end

% 7) If the user has selected 2 markers, then only the pitch angle between 
%    the normal and the vector defined between them can be computed.
if length(segment_select) == 2
    
    % 7.1) Ask the user to select the marker which is placed on the lower 
    %      part of the segment.
    pos_lower_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the lower marker','ListSize',[240 150],...
        'SelectionMode','multiple');

    % 7.2) Extract coordinates from point 1 (lower marker).
    x_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',1,:)')));
    y_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',2,:)')));
    z_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',3,:)')));
    
    
    % 7.3) Extract coordinates from point 2 (upper marker).
    aux = fliplr(segment_select);
    x_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(aux(pos_lower_select)),',1,:)')));
    y_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(aux(pos_lower_select)),',2,:)')));
    z_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(aux(pos_lower_select)),',3,:)')));
    
    % 7.4) Computation of the pitch angle.
    pitch_qs = zeros(1,length(x_coord_1));
    showRT_Plot = 'no'; % Flag which controls if the figure showing the 
                        % real time motion of the two markers is shown.
    showPlot = 'yes';   % Flag which controls if the figure showing the 
                        % computed pitch angle is shown.
    figure
    for i = 1:length(x_coord_1)
        
        % For every time instant, we build a normal vector parallel to the
        % XZ plane which has its origin in PB (lower marker) and finishes
        % in the Z coordinate of PA (upper marker). This defines a right
        % triangle (PA,PB,PC) from which we compute both catheti which are
        % used to compute the angle between vectors PAPB and PBPC. SInce
        % this angle is contained in the XZ plane (i.e. it defines a 
        % rotation around axis Y), so it is defined as the pitch angle.
        
        % We first build the catheti.
        side_a = sqrt((x_coord_2(i) - x_coord_1(i))^2 + (z_coord_2(i) - ...
            z_coord_2(i))^2);
        side_b = sqrt((x_coord_1(i) - x_coord_1(i))^2 + (z_coord_1(i) - ...
            z_coord_2(i))^2);
        
        % And then we apply the arctangent to find the pitch angle.
        pitch_qs(i) = atan2d(side_a,side_b);
        
        % We finally correct the sign of the angle (if needed).
        if x_coord_2(i) < x_coord_1(i)
            pitch_qs(i) = -pitch_qs(i);
        end
        
        % The process is shown graphically in real time (if the flag is
        % activated).
        if strcmpi(showRT_Plot,'yes')
            subplot(2,1,1)
            hold on
            % Lower marker (P1)
            p1 = plot(x_coord_1(i),z_coord_1(i),'.blue','MarkerSize',10);
            % Upper marker (P2)
            p2 = plot(x_coord_2(i),z_coord_2(i),'.red','MarkerSize',10);
            % Point 3 in (x1,y2)
            p3 = plot(x_coord_1(i),z_coord_2(i),'.black','MarkerSize',10);
            % Line joining P1 and P2
            l1 = plot([x_coord_1(i) x_coord_2(i)],[z_coord_1(i) ...
                z_coord_2(i)]);
            % Line joining P1 and P3
            l2 = plot([x_coord_1(i) x_coord_1(i)],[z_coord_1(i) ...
                z_coord_2(i)],'--black');
            pause(1/f)
            delete(p1,p2,p3,l1,l2)
            subplot(2,1,2)
            hold on
            plot(time_qs(i),pitch_qs(i),'MarkerSize',5)
            axis([0 time_qs(end) -180 180])
        end
    end
    
    % 7.5) Correct the sign of the pitch angle to match the GaitWatch 
    %      reference system.
    pitch_qs = -pitch_qs;
    
    % 7.6) In some ocasions the Qualisys data may contain NaN values 
    %      (i.e. there were instans in which all the markers were ocluded). 
    %      Therefore, we need to remove them and interpolate the signal.
    if sum(isnan(pitch_qs)) > 0
        nan_index = isnan(pitch_qs);
        pitch_qs2 = interp1(time_qs(~nan_index),pitch_qs(~nan_index),...
            time_qs(nan_index),'spline');
        pitch_qs(nan_index) = pitch_qs2;
        figure
        plot(time_qs,pitch_qs)
    end
    
end  

% 8) If the user has selected three markers, then both the pitch and roll
%    angles can be computed.
if length(segment_select) == 3    
    
    % 8.1) The user is prompted to select the marker which was placed on 
    %      the upper part of the segment.
    pos_upper_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the upper marker','ListSize',[240 150],...
        'SelectionMode','multiple');
    
    % 8.2) The user is prompted to select the marker which was placed on 
    %      the lower part of the segment.
    pos_lower_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the lower marker','ListSize',[240 150],...
        'SelectionMode','multiple');
    
    % 8.3) The user is prompted to select the marker which was placed on 
    %      the lateral part of the segment.
    pos_side_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the side marker','ListSize',[240 150],...
        'SelectionMode','multiple');
    
    % 8.4) Extract coordinates from point 1 (side marker).
    x_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_select)),',1,:)')));
    y_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_select)),',2,:)')));
    z_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_select)),',3,:)')));

    % 8.5) Extract coordinates from point 2 (upper marker).
    x_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',1,:)')));
    y_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',2,:)')));
    z_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',3,:)')));
    
    % 8.6) Extract coordinates from point 3 (lower marker).
    x_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',1,:)')));
    y_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',2,:)')));
    z_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',3,:)')));
    
    % 8.7) Computation of pitch and roll angles.
    pitch_qs = zeros(1,length(x_coord_1));
    roll_qs = zeros(1,length(x_coord_1));
    showRT_Plot = 'no'; % Flag which controls if the figure showing the 
                        % real time motion of the two markers is shown.
    showPlot = 'yes';   % Flag which controls if the figure showing the 
                        % computed pitch angle is shown.
    figure
    for i = 1:length(x_coord_1)
        
        % For every time instant, we build a normal vector parallel to the
        % XZ plane which has its origin in PB (lower marker) and finishes
        % in the Z coordinate of PA (upper marker). This defines a right
        % triangle (PA,PB,PC) from which we compute both catheti which are
        % used to compute the angle between vectors PAPB and PBPC. Since
        % this angle is contained in the XZ plane (i.e. it defines a 
        % rotation around axis Y), it is defined as the pitch angle. For
        % the roll angle we build a normal vector parallel to the YZ plane
        % which has its origin in P1 (side marker) and finishes in the Z
        % coordinate of PA (upper marker). This defines a right triangle
        % (P1,PA,PD) from which we compute both catheti which are used to
        % compute the angles between vectors P1PA and PAPD. Since this
        % angle is contained in the YZ plane, it is defined as the roll
        % angle.
        
        % We first build the catheti of the pitch triangle.
        side_a = sqrt((x_coord_3(i) - x_coord_2(i))^2 + (z_coord_2(i) - ...
            z_coord_2(i))^2);
        side_b = sqrt((x_coord_3(i) - x_coord_3(i))^2 + (z_coord_2(i) - ...
            z_coord_3(i))^2);
        
        % Now we build the catheti of the roll triangle.
        side_c = sqrt((y_coord_1(i) - y_coord_1(i))^2 + (z_coord_1(i) - ...
            z_coord_2(i))^2);
        side_d = sqrt((y_coord_1(i) - y_coord_2(i))^2 + (z_coord_2(i) - ...
            z_coord_2(i))^2);
        
        % And then we apply the arctangent to find the pitch angle.
        pitch_qs(i) = atand(side_a/side_b);
        roll_qs(i) = atand(side_c/side_d);
        
        % We finally correct the sign of the angle (if needed).
        if x_coord_2(i) < x_coord_3(i)
            pitch_qs(i) = -pitch_qs(i);
        end
        if z_coord_2(i) < z_coord_1(i)
            roll_qs(i) = -roll_qs(i);
        end
        
        % The process is shown graphically in real time (if the flag is
        % activated).
        if strcmpi(showRT_Plot,'yes')
            subplot(4,1,1)
            hold on
            % Lower marker (P3)
            p3 = plot(x_coord_3(i),z_coord_3(i),'.blue','MarkerSize',10);
            % Upper marker (P2)
            p2 = plot(x_coord_2(i),z_coord_2(i),'.red','MarkerSize',10);
            % Point C in (x3,z2)
            pc = plot(x_coord_3(i),z_coord_2(i),'.black','MarkerSize',10);
            % Line joining P3 and P2
            l1 = plot([x_coord_3(i) x_coord_2(i)],[z_coord_3(i) ...
                z_coord_2(i)]);
            % Line joining P3 and PC
            l2 = plot([x_coord_3(i) x_coord_3(i)],[z_coord_3(i) ...
                z_coord_2(i)],'--black');
            axis([1100 1700 0 600])
            pause(1/f)
            delete(p3,p2,pc,l1,l2)
            subplot(4,1,2)
            hold on
            plot(time_qs(i),pitch_qs(i),'MarkerSize',5)
            axis([0 time_qs(end) -100 100])
            subplot(4,1,3)
            hold on
            % side marker (P1)
            p1 = plot(y_coord_1(i),z_coord_1(i),'.green','MarkerSize',10);
            % Upper marker (P2)
            p2 = plot(y_coord_2(i),z_coord_2(i),'.red','MarkerSize',10);
            % Point D in (x1,z2)
            pd = plot(y_coord_1(i),z_coord_2(i),'.cyan','MarkerSize',10);
             % Line joining P1 and P2
            l1 = plot([y_coord_1(i) y_coord_2(i)],[z_coord_1(i) ...
                z_coord_2(i)]);
            % Line joining P1 and PD
            l2 = plot([y_coord_1(i) y_coord_1(i)],[z_coord_1(i) ...
                z_coord_2(i)],'--black');
            axis([1300 1600 300 500])
            pause(1/f)
            delete(p1,p2,pd,l1,l2)
            subplot(4,1,4)
            hold on
            plot(time_qs(i),roll_qs(i),'MarkerSize',5)
            axis([0 time_qs(end) -50 70])
        end
    end

    % 8.7) In some ocasions the Qualisys data may contain NaN values (i.e.
    %      there were instans in which all the markers were ocluded). 
    %      Therefore, we need to remove them and interpolate the signal.
    if sum(isnan(pitch_qs)) > 0
        nan_index = isnan(pitch_qs);
        pitch_qs2 = interp1(time_qs(~nan_index),pitch_qs(~nan_index),...
            time_qs(nan_index),'spline');
        pitch_qs(nan_index) = pitch_qs2;
        figure
        plot(time_qs,pitch_qs)
    end
    
    if sum(isnan(roll_qs)) > 0
        nan_index = isnan(roll_qs);
        roll_qs2 = interp1(time_qs(~nan_index),roll_qs(~nan_index),...
            time_qs(nan_index),'spline');
        roll_qs(nan_index) = roll_qs2;
        figure
        plot(time_qs,roll_qs)
    end
    
    % 8.8) Plot the computed pitch and roll angles.
    if strcmpi(showPlot,'yes')
        figure
        subplot(2,1,1)
        plot(time_qs,pitch_qs)
        xlabel('Time (s)');
        ylabel('Pitch (deg)')
        subplot(2,1,2)
        plot(time_qs,roll_qs)
        xlabel('Time (s)');
        ylabel('Roll (deg)')
    end
end

% 9) If the user has selected four markers, then all the Euler Angles
% (roll, pitch and yaw) can be computed.
if length(segment_select) == 4

    % 9.1) The user is prompted to select the marker which was placed on 
    %      the upper part of the segment.
    pos_upper_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the upper marker','ListSize',[240 150],...
        'SelectionMode','multiple');

    % 9.2) The user is prompted to select the marker which was placed on 
    %      the lower part of the segment.
    pos_lower_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the lower marker','ListSize',[240 150],...
        'SelectionMode','multiple');

    % 9.3) The user is prompted to select the marker which was placed on 
    %      the lateral part of the segment.
    pos_side_a_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the side A marker','ListSize',[240 150],...
        'SelectionMode','multiple');
    
    % 9.4) The user is prompted to select the marker which was placed on 
    %      the lateral part of the segment.
    pos_side_b_select = listdlg('ListString',marker_labels(segment_select),...
        'Name','Select the side B marker','ListSize',[240 150],...
        'SelectionMode','multiple');
    
    % 9.5) Extract coordinates from point 1 (upper).
    x_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',1,:)')));
    y_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',2,:)')));
    z_coord_2 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_upper_select)),',3,:)')));
    
    % 9.6) Extract coordinates from point 2 (lower).
    x_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',1,:)')));
    y_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',2,:)')));
    z_coord_3 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_lower_select)),',3,:)')));
    
    % 9.7) Extract coordinates from point 3 (side A).
    x_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_a_select)),',1,:)')));
    y_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_a_select)),',2,:)')));
    z_coord_1 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_a_select)),',3,:)')));
    
    % 9.8) Extract coordinates from point 3 (side B).
    x_coord_4 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_b_select)),',1,:)')));
    y_coord_4 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_b_select)),',2,:)')));
    z_coord_4 = squeeze(eval(strcat('structName.',names{1},...
        '.Trajectories.Labeled.Data(',...
        num2str(segment_select(pos_side_b_select)),',3,:)')));
    
    % 9.9) Computation of pitch and roll angles.
    pitch_qs = zeros(1,length(x_coord_1));
    roll_qs = zeros(1,length(x_coord_1));
    yaw_qs = zeros(1,length(x_coord_1));
    showRT_Plot = 'no'; % Flag which controls if the figure showing the 
                        % real time motion of the two markers is shown.
    showPlot = 'yes';   % Flag which controls if the figure showing the 
                        % computed pitch angle is shown.
    figure
    for i = 1:length(x_coord_1)
        
        % For every time instant, we build a normal vector parallel to the
        % XZ plane which has its origin in PB (lower marker) and finishes
        % in the Z coordinate of PA (upper marker). This defines a right
        % triangle (PA,PB,PC) from which we compute both catheti which are
        % used to compute the angle between vectors PAPB and PBPC. Since
        % this angle is contained in the XZ plane (i.e. it defines a 
        % rotation around axis Y), it is defined as the pitch angle. For
        % the roll angle we build a normal vector parallel to the YZ plane
        % which has its origin in P1 (side marker) and finishes in the Z
        % coordinate of PA (upper marker). This defines a right triangle
        % (P1,PA,PD) from which we compute both catheti which are used to
        % compute the angles between vectors P1PA and PAPD. Since this
        % angle is contained in the YZ plane, it is defined as the roll
        % angle.
        
        % We first build the catheti of the pitch triangle.
        side_a = sqrt((x_coord_3(i) - x_coord_2(i))^2 + (z_coord_2(i) - ...
            z_coord_2(i))^2);
        side_b = sqrt((x_coord_3(i) - x_coord_3(i))^2 + (z_coord_2(i) - ...
            z_coord_3(i))^2);
        
        % Now we build the catheti of the roll triangle.
        side_c = sqrt((y_coord_1(i) - y_coord_1(i))^2 + (z_coord_1(i) - ...
            z_coord_2(i))^2);
        side_d = sqrt((y_coord_1(i) - y_coord_2(i))^2 + (z_coord_2(i) - ...
            z_coord_2(i))^2);
        
        % And then we apply the arctangent to find the pitch angle.
        pitch_qs(i) = atand(side_a/side_b);
        roll_qs(i) = atand(side_c/side_d);
        
        % We finally correct the sign of the angle (if needed).
        if x_coord_2(i) < x_coord_3(i)
            pitch_qs(i) = -pitch_qs(i);
        end
        if z_coord_2(i) < z_coord_1(i)
            roll_qs(i) = -roll_qs(i);
        end
        
        % The process is shown graphically in real time (if the flag is
        % activated).
        if strcmpi(showRT_Plot,'yes')
            subplot(4,1,1)
            hold on
            % Lower marker (P3)
            p3 = plot(x_coord_3(i),z_coord_3(i),'.blue','MarkerSize',10);
            % Upper marker (P2)
            p2 = plot(x_coord_2(i),z_coord_2(i),'.red','MarkerSize',10);
            % Point C in (x3,z2)
            pc = plot(x_coord_3(i),z_coord_2(i),'.black','MarkerSize',10);
            % Line joining P3 and P2
            l1 = plot([x_coord_3(i) x_coord_2(i)],[z_coord_3(i) ...
                z_coord_2(i)]);
            % Line joining P3 and PC
            l2 = plot([x_coord_3(i) x_coord_3(i)],[z_coord_3(i) ...
                z_coord_2(i)],'--black');
            axis([1100 1700 0 600])
            pause(1/f)
            delete(p3,p2,pc,l1,l2)
            subplot(4,1,2)
            hold on
            plot(time_qs(i),pitch_qs(i),'MarkerSize',5)
            axis([0 time_qs(end) -100 100])
            subplot(4,1,3)
            hold on
            % side marker (P1)
            p1 = plot(y_coord_1(i),z_coord_1(i),'.green','MarkerSize',10);
            % Upper marker (P2)
            p2 = plot(y_coord_2(i),z_coord_2(i),'.red','MarkerSize',10);
            % Point D in (x1,z2)
            pd = plot(y_coord_1(i),z_coord_2(i),'.cyan','MarkerSize',10);
             % Line joining P1 and P2
            l1 = plot([y_coord_1(i) y_coord_2(i)],[z_coord_1(i) ...
                z_coord_2(i)]);
            % Line joining P1 and PD
            l2 = plot([y_coord_1(i) y_coord_1(i)],[z_coord_1(i) ...
                z_coord_2(i)],'--black');
            axis([1300 1600 300 500])
            pause(1/f)
            delete(p1,p2,pd,l1,l2)
            subplot(4,1,4)
            hold on
            plot(time_qs(i),roll_qs(i),'MarkerSize',5)
            axis([0 time_qs(end) -50 70])
        end
    end
    
    % 9.10) Correct the sign of the pitch angle to match the GaitWatch 
    %      reference system.
    pitch_qs = -pitch_qs;
    % 9.11) Plot the computed angles.
    if strcmpi(showPlot,'yes')
        figure
        subplot(2,1,1)
        plot(time_qs,pitch_qs)
        xlabel('Time (s)');
        ylabel('Pitch (deg)')
        subplot(2,1,2)
        plot(time_qs,roll_qs)
        xlabel('Time (s)');
        ylabel('Roll (deg)')
    end
    
end

