function [header,first_line,abta,count,col_first_block] = read_laufband_header(filename)
% das Programm liest den header von textfiles, die von zebris Ganganalyse stammen
% header{1}: Name der VP
% header{2}: Application
% header{3}: data_labels
% header{4}: Calibration Zebris
% header{5}: Creation date (string)
% abta : % Abtastintervall in ms
% count : Anzahl von Messungen/datenblöcken
% first_line: Beginn des datenblocks im textfile
% calib_zebris: Die  Kalibrierwerte sind die x,y,z Koordinaten der Marker
% col_first_block= no of col first block
% während der Kalibrierung – geben demnach die Neutralposition des Probanden an.
%
% Die Werte, welche unter „Verschieben,cell“ angeben sind, stellt die Positionsänderung des Laufbandgurts dar.

 

fid=fopen(filename);
disp(filename)
C = textscan(fid, '%s',19,'Delimiter','\n','whitespace', '');   % header lesen 19 Zeilen


ind=find(strncmp('Per',C{1},3));                          % Name
if ~isempty(ind)
        nn=textscan(char(C{1}(ind)),'%s');
        nn=char(nn{:}); nn(1,:) = []; [r,c] = size(nn);                    % contains 'Person'
        header{1}=reshape(nn',1,r*c);
end

ind=find(strncmp('App',C{1},3));                          % Application: typ des Versuches und Typ des Datenfiles !!
nn=textscan(char(C{1}(ind)),'%s%s');
header{2}=char(nn{2});

ind=find(strncmp('Cre',C{1},3));                          % Creation date
nn=textscan(char(C{1}(ind)),'%s%s%s%s');
% disp(char(C{1}(ind)))
header{5}=[char(nn{3}),' ',char(nn{4})];
% disp(header{5})

ind=find(strncmp('Fre',C{1},3));                          % Abtastintervall in ms
nn=textscan(char(C{1}(ind)),'%s%n');
abta=1000/nn{2}; 

ind=find(strncmp('Cou',C{1},3));                           % zeile 15: count, Anzahl von Messungen/Datenblöcken
nn=textscan(char(C{1}(ind)),'%s%n');
count=nn{2};

ind=find(strncmp('Tim',C{1},3));                           % zeile 17: Data labels bzw Anzahl der Spalten im ersten block bei application trplatf
if strcmp(header{2}, 'trplatf') || strcmp(header{2}, 'statplat')
     data_labels=textscan(char(C{1}(ind)),'%s','delimiter','\t');
     col_first_block=length(data_labels{1});                                  % Spalten im ersten block
     header(3)=data_labels;
% elseif strcmp(header{2},'statplat')
%      col_first_block=1;                                               % dabei nur eine Spalte
end

ind=find(strncmp('Cal',C{1},3));                           % calibration for Zebris 18 Werte !! cal ist die vorletzte zeile bevor daten kommen
if strcmp(header{2}, 'trplatf')
     calib=textscan(char(C{1}(ind)),'%n','delimiter','\t','treatAsEmpty','Calibr');  % calibrationsdaten zebris Armswing-Marker
     header{4}=calib{1}(~isnan(calib{1}));
else header{4}=0;
end
first_line=ind+1;                                          % data begin in first line
fclose(fid);
