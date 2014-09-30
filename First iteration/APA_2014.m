% APA_2014.m
% basiert auf: gehen_laufband.m
% program reads textfiles made by Zebris pressure-sensitive footplate
% For evaluation of APA data obtained in Munich 2014

clear
[bb,aa]=butter(3,0.3);
% Maybe other simplier form to do : [data_file,data_path] = uigetfile('*.txt', 'Select *.txt - file');
data_path='C:\A-DATA\';

% Data file selection with a dialog box (only the text file)
[data_file,data_path] = uigetfile([data_path,'*.txt'], 'Select *.txt - file');


%Reading of data_file
[header,first_line,abta,count,col_first_block] = read_laufband_header([data_path,'\',data_file]);
 % header: file name, application, data labels, calibration and creation
 % date.
 % first_line:Beginning of the data block in the text file.
 % abta: sampling time (ms).
 % count:Number of measurements/ data .
 % col_first_block: no of col first block.


% ********** Read first block of data ********************************************


% First block contains 5 columns:
%     'Time,ms'
%     'Li Vorfuﬂ,N'
%     'Li R¸ckfuﬂ,N'
%     'Re Vorfuﬂ,N'
%     'Re R¸ckfuﬂ,N'

% Read formatted data from text file or string
% "fid" is a file identifier that we obtain with "fopen".
    
fid=fopen([data_path,'\',data_file]);

% Specification of the format of the data fields (double, %n) an cycles. It
% skips the firsts lines of the data, and then reads the remaining data.
data = textscan(fid,'%n',col_first_block*count,'headerlines',first_line-1);

% Reshape the data.
data=reshape(data{1},col_first_block,count)';

an_en=any(data(:,2:5),2);                                  % see if there are data available 
an_en=find(diff(an_en) ~= 0)+1;                            % where do the data begin and end ??

% Show the data
figure(1)
plot(data(an_en(1):an_en(2),1),data(an_en(1):an_en(2),2:5))
legend(header{3}(2:5))
xlabel(header{3}(1))
ylabel('N')



% % read some more details about the following blocks
% 
C2 = textscan(fid, '%s',10,'Delimiter','\n','whitespace', '');   % read 10 lines which contain information concerning the following data blocks
nn=textscan(char(C2{1}(5)),'%s%n');                              % find out no of colums per block
colums=nn{2};
nn=textscan(char(C2{1}(6)),'%s%n');                             % find out no of lines per block 
lines=nn{2};
data2=zeros(lines,colums);
format=['%s',repmat('%n',1,colums)];                            % format for reading data, first value is a string



% read data block by block and evaluate only those which contain data
% only to define the midline 
m=zeros(colums,1);           

for b=1:an_en(2)
     CX = textscan(fid, '%s',8,'Delimiter','\n','whitespace', '');   % read 8 lines which contain some irrelevant information
     C3 = textscan(fid,format,lines);
     if b >= an_en(1)
        for c=1:colums
            m(c) = m(c) + mean(C3{c+1});                             % sum up the mean of all columns
        end
     end
end

fclose(fid);

% define midline between both feet

[pks,locs] = findpeaks(m,'minpeakdistance',10,'minpeakheight',100);
if length(pks) == 2
    [~,ml] = min(m(locs(1):locs(2)));
    ml = ml + locs(1) - 1;
else
    disp('midline could not be defined')
    beep
end

figure(2)
plot(m)
title('Midline between both feet')
xlabel('pressure in r-l-dimension')
text(locs, pks,'+','VerticalAlignment','Bottom','HorizontalAlignment','Center','FontSize',16,'col','r')
text(ml,max(m)/10,'\downarrow','VerticalAlignment','Bottom','HorizontalAlignment','Center','FontSize',16,'col','r')





% read data one more time, plot them and evaluate

n=an_en(2)-an_en(1)+1;
data3=zeros(n,3,3);   % 3 values per frame(x,y,N), for 3 feet: No 1 is right foot 2 is left foot 3 is both feet
% Werte: 1: cog x (right-left)
%        2: cog y  (ant-post)
%        3: posterrior margin (only y)
%        4: anteriour margin, der belastet ist (nur y)
%        5: overall pressure 

%              prepare graphic
data2 = zeros(lines,colums);
figure(1),clf
axes('pos',[ 0.1300    0.1100    0.7750    0.8150]);
ph=pcolor(data2);
axis equal
axis([1 colums 1 lines])
drawnow
th=text(60,45,num2str(b));

lh=line([ml ml],[1,lines],'color',[1 1 1]);

cogr_h=line([20 25] ,[20 20] ,'color',[1 1 1],'linewidth',2);                   % graphic marker COG right
cogl_h=line([20 25] ,[20 20] ,'color',[1 1 1],'linewidth',2);
cogb_h=line([20 25] ,[20 20] ,'color',[1 0 0],'linewidth',3);                   % graphic marker COG both
graph_handles=[cogr_h cogl_h cogb_h];


           
           
fid=fopen([data_path,'\',data_file]);
dummy_data = textscan(fid,'%n',col_first_block*count,'headerlines',first_line-1);
C2 = textscan(fid, '%s',10,'Delimiter','\n','whitespace', '');


for b=1:an_en(1)-1                                                   % skip these data they contain no values
     CX = textscan(fid, '%s',8,'Delimiter','\n','whitespace', '');   % read 8 lines which contain some irrelevant information
     C3 = textscan(fid,format,lines);
end


for b=an_en(1):an_en(2)
    c=b-an_en(1)+1;                                                 % index for data3
    CX = textscan(fid, '%s',8,'Delimiter','\n','whitespace', '');   % read 8 lines which contain some irrelevant information
    C3 = textscan(fid,format,lines);
    for cc=1:colums
       data2(:,cc)=C3{cc+1};
    end
    set(ph,'cdata',data2)
    set(th,'string',num2str(b))
    for side=1:3                                                              % Seite 1 rechts, Seite 2 links seite 3: both
           switch side
               case 1
                    an=1; en=ml;
               case 2
                    an=ml+1; en=colums;                                            % colums = x=(64); lines = y=(40); size(data2) = 40    64
               case 3
                   an=1; en=colums;
           end

           m=sum(data2(:,an:en));                                           % Colum vektor! Berechnung der x-Koordinate des Schwerpunktes
           if any(m)                                                           % pr¸fen, ob ¸berhaupt in diesem Feld Daten>0 sind 
               data3(c,1,side)= sum(m .* (an:en)) ./ sum(m);                  % Ergebnisarray belegen x-Koordinate des Schwerpunktes
               m=sum(data2(:,an:en),2);                                          % Berechnung der y-Koordinate des Schwerpunktes
               data3(c,2,side)= sum(m .* (1:lines)') ./ sum(m);
               data3(c,3,side)=sum(m);                                                  % hintersten Fuﬂdruckpunkt (Fersenposition) berechnen (nur y)
               set(graph_handles(side),'xdata',[data3(c,1,side)-2 data3(c,1,side)+2],'ydata',[data3(c,2,side) data3(c,2,side)],'vis','on')
           else
                 set(graph_handles(side),'vis','off')                         % auf dieser Seite kein Fuﬂ auf dem Boden
           end
            
    end
    drawnow
end
   
fclose(fid);

data3(data3==0) = NaN;

l={'r foot','l foot','b feet'};

figure(3)
subplot(3,1,1)
plot(data(an_en(1):an_en(2),1), squeeze(data3(:,1,:)),'linewidth',1.5)
legend(l)
title('lateral COG excursions')
xlabel('Time [ ms ]')
ylabel('< r     l >')

subplot(3,1,2)
plot(data(an_en(1):an_en(2),1), squeeze(data3(:,2,:)),'linewidth',1.5)
legend(l)
title('ant-post COG excursions')
xlabel('Time [ ms ]')
ylabel('< p     a >')

subplot(3,1,3)
plot(data(an_en(1):an_en(2),1), squeeze(data3(:,3,:)),'linewidth',1.5)
legend(l)
title('Force')
xlabel('Time [ ms ]')
ylabel('[ N ]')




axes('position', [0 0 1 1],'visible', 'off')
text(0.5,0.05,strrep(data_file,'_',' '),'units','normalized','horizontalAlignment','center','VerticalAlignment','top','Fontsize',12)







