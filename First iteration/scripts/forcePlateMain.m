% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ------------------- ZEBRIS FORCE PLATE DATA ANALYSIS --------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai B�tzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Ver�nica Torres (2): 
%                   |_ vts24@correo.ugr.es 
%                 - Dr. Eng. Alberto Olivares (3): 
%                   |_ aolivares@ugr.es
%                 - Robin Weiss (4): 
%                   |_ mail@robinweiss.de
%
% * Affiliation: (1) Laboratory of Motion Analysis and Deep Brain 
%                    Stimulation, Department of Neurology, Klinikum 
%                    Grosshadern, Munich, Germany.
%                (2) Master in Telecommunication Engineering, University of 
%                    Granada, Granada, Spain, (student).
%                (3) Signal Processing and Biomedical Applications group,
%                    Department of Signal Theory, Telematics and
%                    Communications, University of Granada, Granada, Spain.
%                (4) Bachelor in Electrical Engineering, University of 
%                    Applied Sciences of Munster, Munster, Germany, 
%                    (student).
%
% * Last modification: 08/10/2014
% -------------------------------------------------------------------------
% INFORMATION: This file contains the routine to read and plot the data
% coming from a Zebris FDM-S Multifunction Force-measuring Plate used to
% measure gait initiation and balance of Parkinson patients. This is a
% modified version of the 'gehen_laufband.m' file implemented by (?? ask 
% this).
% -------------------------------------------------------------------------

% Delete workspace.
clear all

% Build Butterworth filter.
[bb, aa] = butter(3, 0.3);

% Select data file with a dialog box (only .txt files).
[filename, filepath] = uigetfile('*.txt', ...
    'Select Force Plate data file (.txt)', '../data/ForcePlate');

% Read data file.
[header, first_line, abta, count, col_first_block] = ...
    read_fp_header(fullfile(filepath,filename));

% Read first block of data. The first block contains 5 columns:
%   (1) Time (ms).
%   (2) 'Li Vorfu�,N':  Front left foot pressure sensor.
%   (3) 'Li R�ckfu�,N': Back left foot pressure sensor.
%   (4) 'Re Vorfu�,N':  Front right foot pressure sensor.
%   (5) 'Re R�ckfu�,N': Back right foot pressure sensor.

% Read formatted data from text file. As always "fid" is a file identifier 
% that we obtain with "fopen".
fid = fopen(fullfile(filepath, filename));

% Specification of the format of the data fields (double, %n) an cycles. It
% skips the firsts lines of the data, and then reads the remaining data.
data = textscan(fid, '%n', col_first_block * count, 'headerlines', ...
    first_line - 1);

% Reshape the data. This converts the data from a single column block to a
% five column block.
data = reshape(data{1}, col_first_block, count)';

% Find the starting and ending point of data. This removes all the zero
% readings, i.e. the time instants in which the patient is not walking on 
% the force plate.
start_end = any(data(:, 2:5), 2);                                  
start_end = find(diff(start_end) ~= 0)+1; 
data_start = start_end(1);
data_end = start_end(2);

% Plot the force data.
figure(1)
plot(data(data_start:data_end, 1), data(data_start:data_end, 2:5))
legend(header{3}(2:5))
xlabel(header{3}(1))
ylabel('Force (N)')

% -------------------------------------------------------------------------
% Read some more details about the following blocks.
% -------------------------------------------------------------------------
% Read 10 lines which contain information concerning the following data 
% blocks.
C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');   

% Find out no of columns per block.
nn = textscan(char(C2{1}(5)), '%s%n');                              
columns = nn{2};

% Find out no of lines per block 
nn = textscan(char(C2{1}(6)), '%s%n');                             
lines = nn{2};

% Set format to read the data, first value is a 'string' variable, the rest
% are 'double' variables.
format = ['%s', repmat('%n', 1, columns)];                            

% Read data block by block and evaluate only those which contain data
% only to define the midline. 
m = zeros(columns,1);           

for b = 1:data_end
     % Read 8 lines which contain some irrelevant information.
     CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
     C3 = textscan(fid, format, lines);
     if b >= data_start
        for c = 1:columns
            % Sum up the mean of all columns.
            m(c) = m(c) + mean(C3{c + 1});                             
        end
     end
end

% Close data file.
fclose(fid);

% -------------------------------------------------------------------------
% Definition of the midline between both feet.
% -------------------------------------------------------------------------
[pks, locs] = findpeaks(m, 'minpeakdistance', 10, 'minpeakheight', 100);
if length(pks) == 2
    [~,ml] = min(m(locs(1):locs(2)));
    ml = ml + locs(1) - 1;
else
    disp('midline could not be defined');
    beep;
end

figure(2);
plot(m);
title('Midline between both feet');
xlabel('pressure in r-l-dimension');
text(locs, pks, '+', 'VerticalAlignment', 'Bottom', ...
    'HorizontalAlignment', 'Center', 'FontSize', 16, 'col', 'r');
text(ml, max(m) / 10, '\downarrow', 'VerticalAlignment', 'Bottom', ...
    'HorizontalAlignment','Center','FontSize',16,'col','r')

% -------------------------------------------------------------------------
% Read data one more time, plot them and evaluate.
% -------------------------------------------------------------------------

n = data_end - data_start + 1;

% 3 values per frame(x,y,N), for 3 feet: No 1 is right foot, 2 is left foot
% and 3 is both feet.
data3 = zeros(n,3,3); 

% Values: (1): COG (Center of Gravity) X (right-left).
%         (2): COG Y (anteroposterior i.e. front-to-back).
%         (3): Posterior margin (only Y).
%         (4): Anterior margin (only Y is loaded).
%         (5): Overall pressure.

% Prepare figure.
data2 = zeros(lines, columns);
figure(1);
clf;
axes('pos', [ 0.1300    0.1100    0.7750    0.8150]);
ph = pcolor(data2);
axis equal;
axis([1 columns 1 lines]);
drawnow;
th = text(60,45,num2str(b));

% Draw central line
lh = line([ml ml], [1, lines], 'color',[1 1 1]);

% Initialize graphic marker COG (right).
cogr_h = line([20 25], [20 20], 'color', [1 1 1], 'linewidth', 2);

% Initialize graphic marker COG (left).
cogl_h = line([20 25] ,[20 20] ,'color',[1 1 1],'linewidth',2);

% Initialize graphic marker COG (both).
cogb_h = line([20 25] ,[20 20] ,'color',[1 0 0],'linewidth',3);                   
graph_handles = [cogr_h cogl_h cogb_h];

           
fid = fopen(fullfile(filepath,filename));
dummy_data = textscan(fid, '%n', col_first_block * count, 'headerlines',...
    first_line - 1);
C2 = textscan(fid, '%s', 10, 'Delimiter', '\n', 'whitespace', '');

% Skip these data they contain no values
for b = 1:data_start - 1  
     % Read 8 lines which contain some irrelevant information
     CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
     C3 = textscan(fid, format, lines);
end


for b = data_start:data_end
    % index for data3
    c = b - data_start + 1; 
    % read 8 lines which contain some irrelevant information
    CX = textscan(fid, '%s', 8, 'Delimiter', '\n', 'whitespace', '');   
    C3 = textscan(fid, format, lines);
    for cc = 1:columns
       data2(:, cc) = C3{cc + 1};
    end
    set(ph, 'cdata', data2)
    set(th, 'string', num2str(b))
    % Side 1 right, side 2 left, side 3 both.
    for side = 1:3                                                              
        switch side
            case 1
                an = 1; 
                en = ml;
            case 2
                an = ml + 1;
                en = columns;                                            
            case 3
                an = 1; 
                en = columns;
        end
        % Colum vektor! Berechnung der x-Koordinate des Schwerpunktes
        m = sum(data2(:, an:en));   
        
        % pr�fen, ob �berhaupt in diesem Feld Daten>0 sind 
        if any(m)    
            % Ergebnisarray belegen x-Koordinate des Schwerpunktes
            data3(c,1,side) = sum(m .* (an:en)) ./ sum(m); 

            % Berechnung der y-Koordinate des Schwerpunktes
            m = sum(data2(:, an:en), 2);                                          
            data3(c, 2, side) = sum(m .* (1:lines)') ./ sum(m);

            % hintersten Fu�druckpunkt (Fersenposition) berechnen (nur y)
            data3(c, 3, side) = sum(m);                                                  
            set(graph_handles(side), 'xdata', [data3(c, 1, side) - 2 ...
               data3(c, 1, side) + 2], 'ydata', [data3(c, 2, side) ...
               data3(c, 2, side)], 'vis', 'on');
        else
            % auf dieser Seite kein Fu� auf dem Boden
            set(graph_handles(side),'vis','off');                         
        end
            
    end
    drawnow;
end

% Close data file.
fclose(fid);

%--------------------------------------------------------------------------
% Plot all data3.
%--------------------------------------------------------------------------
data3(data3==0) = NaN;

figure(3)
% Set legend
l = {'r foot','l foot','b feet'};

% Plot right-left COG
subplot(3, 1, 1);
plot(data(data_start:data_end, 1), squeeze(data3(:, 1, :)), 'linewidth',...
    1.5);
legend(l);
title('lateral COG excursions');
xlabel('Time [ ms ]');
ylabel('< r     l >');

% Plot front-back COG
subplot(3,1,2);
plot(data(data_start:data_end, 1), squeeze(data3(:, 2, :)), 'linewidth',...
    1.5);
legend(l);
title('ant-post COG excursions');
xlabel('Time [ ms ]');
ylabel('< p     a >');

% Plot force.
subplot(3,1,3);
plot(data(data_start:data_end, 1), squeeze(data3(:, 3, :)), 'linewidth',...
    1.5);
legend(l);
title('Force');
xlabel('Time [ ms ]');
ylabel('[ N ]');

axes('position', [0 0 1 1],'visible', 'off')
% Print filename on the figure.
text(0.5, 0.05, strrep(filename, '_', ' '), 'units', 'normalized', ...
    'horizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'Fontsize',12);

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% END OF FORCEPLATEMAIN FILE
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\