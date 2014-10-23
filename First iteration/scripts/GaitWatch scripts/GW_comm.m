% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% ||||||||||||||| GAITWATCH'S USB COMMUNICATION ROUTINE |||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%
% GW_COMM communicates with GaitWatch over USB. A list of the possible 
% actions is shown in the command line. Current version's possible
% actions are:
%     - 'Read and change setup': It gives the following options:
%           - 1 Fileprefix                     
%           - 2 Trigger intevall               
%           - 3 Trigger duration               
%           - 4 Sampling frequency             
%           - 5 Settings for magnet            
%           - 6 Settings for accelerometer     
%           - 7 Settings for gyros             
%           - 9 ready -> send to GaitWatch     
%           - 10 Exit withouth updating settings
%     - 'Get all the file headers': Shows a list showing the information of 
%       all the data files stored in the memory card. The user is then
%       asked whether he wants to load a data file (or more than one). 
%       If so, the user needs to input the ID of the initial and final
%       files to be loaded. If only one file is to be loaded, initial and
%       final ID should be the same.
%     - 'Get time and date': (Here goes a brief explanation)
%     - 'Set time and date': (Here goes a brief explanation)
%
% 
%
% |_ 'data': Data matrix containing all the measurements of the sensors.
% |_ 'FileHeader': Header of the file containing the following information:
%                 |_ 'length': 
%                 |_ 'version': Version of GaitWatch's firmware.
%                 |_ 'NChannels': Number of data channels (see main.m file
%                    for more information about the channels).
%                 |_ 'time_begin': Start time of data gathering:
%                       - time_begin(1): second.
%                       - time_begin(2): minute.
%                       - time_begin(3): hour.
%                       - time_begin(5): day.
%                       - time_begin(6): month.
%                       - time_begin(7): year.
%                 |_ 'time_end': End time of data gathering:
%                       - time_end(1): second.
%                       - time_end(2): minute.
%                       - time_end(3): hour.
%                       - time_end(5): day.
%                       - time_end(6): month.
%                       - time_end(7): year.
%                 |_ 'Nsamples': Number of samples in the data matrix.
%                 |_ 'firstBlock': address in memory of the first block of
%                    data.
%                 |_ 'lastBlock': address in memory of the last block of
%                    data.
%                 |_ 'FileNumber': ID of the file (different to the number
%                    determining the order in the file list which is shown 
%                    to the user.
%                 |_ 'Fileprefix':
%                 |_ 'TringInt':
%                 |_ 'TrigDur':
%                 |_ 'SampFreq': Sampling frequency.
%                 |_ 'RecMode':
%                 |_ 'Mag_set': configuration parameters of the
%                    magnetometer (see magnetometer's datasheet).
%                 |_ 'ACC_set': configuration parameters of the
%                    acceleremeter (see accelerometer's datasheet).
%                 |_ 'IMU_set': 
%                 |_ 'ChannelNames': A list containing the name of each 
%                 data channel.
%                 |_ 'ScaleFactors': 
%                 |_ 'Offsets':
%
% ************************************************************************* 
% -------------------------------------------------------------------------
% Authors: Kai Bötzel & Alberto Olivares.
% Entity: Ludwig-Maximilians Universität München & Universidad de Granada.
% Last modification: 07/11/2013.
% -------------------------------------------------------------------------


% message from avr to matlab

%    0 1. testbyte, immer 0xAAAA
% 	 1 2. Datentyp =   80 = uint8,
%                     81 = int8, (signed)
%                     160 = uint16,
%                     161 = int16 (signed)
%                     32 = uint32
%                     132 = int32
%                     100 = character (ASCII)
%                      0 = keine Daten, nur Matlab message
%                      20 = File header (char)
%  	2  3. wieviele Werte werden von dem genannnten Datentyp insgesamt gesendet uint32 high byte
% 	3  4. wieviele Werte werden von dem genannnten Datentyp insgesamt gesendet uint32 low byte
% 	4  5. wieviele Werte (nicht byte) hat ein gesendeter Block (kann variieren)
%   5  6. Anzahl an Kanälen
%   6  7. confirmation expected = 0, keine Antwort nötig, 1: AVR wartet auf  Antwort
%   7  8. Kann unterschiedlich belegt werden Dauer eines Zyklus,  Abtastfrequenz o.ä. 
%   nach dieser message kommen die daten

% Matlab arbeitet immer mit Angaben von Werten nicht bytes
% avr arbeitet immer mit bytes

% s = serial('COM16','BaudRate',57600,'Timeout',20);
% set(s,'InputBufferSize',5000)
% fopen(s);

if isempty(instrfind)
    s = serial('COM6','BaudRate',512000,'Timeout',10);
    set(s,'InputBufferSize',4096)
    fopen(s);
end


m='The specified amount of data was not returned within the Timeout period.';
NPlot= 200;           % Länge der Plotdaten
len_mess=8;           % Länge matlab-message in integer
% set(s,'Timeout',30);
OK=1;

txt={' ';...
     'Menue:';...
       '   190    Read & change setup';...
     % '   200    Online Data transfer';...
     % '   150    Accelerometer Monitor';...     
     % '   210    Record to SD-card';...
     '   215    get all file headers';...          % task 220 ist vergeben !!
     % '   125    read error log';...
     % '   240    delete last file';...     
     % '   250    Reset SD card / delete all data';...
       '   30     Get Time&Date';...
       '   31     Set Time&date';...
     % '   120    Check Battery';...
     '   99     Quit';...
     };
 

while OK
   disp(char(txt))
   task=input('Bitte eingeben: ');
   
   count=get(s,'bytesavailable'); 
   while count > 0
       flush(count,s)
       count=get(s,'bytesavailable');        
   end
   
   
   if task ~= 99
      fwrite(s,[hex2dec('AA'),task],'uint8')                                        % send command
   end
   
   
   switch task
       
      case 30                                                           % get clock
          [mat_mess,count,msg] = fread(s, len_mess, 'uint16');
          if isempty(msg)  && mat_mess(1) ==  hex2dec('AAAA')
             [TimeDate,count,msg] = fread(s, 7, 'uint8');
              disp(sprintf('Datum : %02d.%02d.%d       Uhrzeit : %02d:%02d:%02d',...
              TimeDate(5),TimeDate(6),TimeDate(7)+2000,TimeDate(3),TimeDate(2),TimeDate(1)))
          end
          
          % erkennen, wann die uhr einen reset gemacht hat
%           ResetDate=now - datenum([TimeDate(7) TimeDate(6) TimeDate(5) TimeDate(3),TimeDate(2),TimeDate(1)]); %datevector!
%           disp(datestr(ResetDate))
%  
          
       case 31                                                           % set clock
          TimeDate=uint8(zeros(1,7));
          temp=round(datevec(now));
          TimeDate(1:3)=temp(6:-1:4);
          TimeDate(5:6)=temp(3:-1:2);
          TimeDate(7)=temp(1)-2000;
          fwrite(s,TimeDate,'uint8')
          disp('data were sent!')
          [mat_mess,count,msg] = fread(s, len_mess, 'uint16');           % Matlab wartet auf Antwort (Textstring)
          [mess,count,msg] = fread(s,mat_mess(5),'char');
          disp(char(mess'))
          if ~isempty(msg), disp(msg),end
          
       
       case 215                                                                     % 215: read file directory
           [mat_mess,count,msg] = fread(s, len_mess, 'uint16');                     % Mat mess einmal lesen
            fwrite(s,[hex2dec('AA'),100],'uint8')                % confirm
            
            SD_ADDR = typecast(uint16([mat_mess(3) mat_mess(2)]),'uint32');         % First block of next file on SD card
            SD_Max_Blocks = typecast(uint16([mat_mess(5) mat_mess(4)]),'uint32');   % First block of next file on SD card
            blocksize   = mat_mess(6);                                              % Länge fileheader
            Files_on_SD = mat_mess (7);
            fprintf('\nMaximum Blocks (512 bytes) on SD-card        : %d \n',SD_Max_Blocks)
            fprintf('Address of next Block on SD on SD-card       : %d \n',SD_ADDR)
            fprintf('Space occupied on SD-card                    : %-5.2f %% \n',double(SD_ADDR) / double(SD_Max_Blocks) * 100)
            fprintf('Number of files on SD-card                   : %d \n', Files_on_SD)
            fprintf('Maximum number of files in SD-card directory : %d \n\n\n', mat_mess(8))

           %   ***************            prepare matlab array for retrieval of SD-data    **********************

           SD_adresses = uint32(zeros(Files_on_SD+1,1));


           %   ***************            Read all file headers    **********************
           clear FileHeader
            n_files = 0;
            while n_files < Files_on_SD
                [fh,count,msg] = fread(s, blocksize, 'uint8');
                fwrite(s,[hex2dec('AA'),100],'uint8')                              % confirm
                FileHeader.length = typecast(uint8(fh(1:2)),'uint16');          % file header length in bytes
                FileHeader.version    = fh(3); 
                FileHeader.NChannels  = fh(4);
                FileHeader.time_begin = fh(5:11)';
                FileHeader.time_end = fh(12:18)';
                FileHeader.NSamples = typecast(uint8(fh(19:22)),'uint32');
                FileHeader.firstBlock = typecast(uint8(fh(23:26)),'uint32');
                FileHeader.lastBlock = typecast(uint8(fh(27:30)),'uint32');
                FileHeader.FileNumber = typecast(uint8(fh(31:32)),'uint16');
                FileHeader.Fileprefix = fh(33:42)';
                FileHeader.TrigInt = typecast(uint8(fh(43:44)),'uint16');
                FileHeader.TrigDur = typecast(uint8(fh(45:46)),'uint16');
                
                FileHeader.FileName = [char(FileHeader.Fileprefix),sprintf('%04d',FileHeader.FileNumber)];
                n_files = n_files + 1;
                SD_adresses(n_files) = FileHeader.firstBlock;                % store information of file headers  for retrieval of data from SD card 
                dauer=FileHeader.time_end(1:3)-FileHeader.time_begin(1:3); 
                for n=1:2
                    if dauer(n) < 0
                        dauer(n)=60+dauer(n);
                        dauer(n+1)=dauer(n+1)-1;
                    end
                end
               %   ***************            display information of file headers    **********************                    
                    fprintf('%d.)  %s   Datum : %02d.%02d.%d    Uhrzeit : %02d:%02d:%02d    Dauer : %02d:%02d:%02d \n',...
                        n_files,... 
                        FileHeader.FileName,...
                        FileHeader.time_begin(5),...
                        FileHeader.time_begin(6),...
                        FileHeader.time_begin(7)+2000,FileHeader.time_begin(3),...
                        FileHeader.time_begin(2),FileHeader.time_begin(1),...
                        dauer(3),...
                        dauer(2),...
                        dauer(1));
                

            end 
            SD_adresses(n_files+1) = SD_ADDR;                                         % First block of next file on SD card
                                
            
            while input('Want to get data from SD-Card ? (y/n) : ','s') == 'y'
               first_file= input('First file number : ');
               last_file = input('Last file number  : ');
               blocks_to_read = SD_adresses(last_file+1) - SD_adresses(first_file);       % number of data blocks to read
               blocks_read = 0;
               h = waitbar(0,'Please wait...');
               for n_files = first_file:last_file             
                    fwrite(s,[hex2dec('AA'),220],'uint8')                        % task 220 read data fron SDcard
                    fwrite(s,[SD_adresses(n_files) (SD_adresses(n_files+1)-1)],'uint32')  % send blocks to read
                    % read fileheader (first block)
                    [fh,count,msg] = fread(s, 512, 'uint8');                             % first block contains fileheader
                    blocks_read = blocks_read+1;
                    fwrite(s,[hex2dec('AA'),100],'uint8')                                      % confirm
                    clear FileHeader
                    FileHeader.length = typecast(uint8(fh(1:2)),'uint16');          % file header length in bytes
                    FileHeader.version    = fh(3); 
                    FileHeader.NChannels  = fh(4);
                    FileHeader.time_begin = fh(5:11)';
                    FileHeader.time_end = fh(12:18)';
                    FileHeader.NSamples = typecast(uint8(fh(19:22)),'uint32');
                    FileHeader.firstBlock = typecast(uint8(fh(23:26)),'uint32');
                    FileHeader.lastBlock = typecast(uint8(fh(27:30)),'uint32');
                    FileHeader.FileNumber = typecast(uint8(fh(31:32)),'uint16');
                    FileHeader.Fileprefix = fh(33:42)';
                    FileHeader.Fileprefix(FileHeader.Fileprefix==0)=[];
                    FileHeader.TrigInt = typecast(uint8(fh(43:44)),'uint16');
                    FileHeader.TrigDur = typecast(uint8(fh(45:46)),'uint16');
                    FileHeader.SampFreq = typecast(uint8(fh(47:48)),'uint16');
                    FileHeader.RecMode = fh(49);
                    FileHeader.Mag_set = fh(50:53)';
                    FileHeader.ACC_set = fh(54:57)';
                    FileHeader.IMU_set = fh(58:61)';
                    
                    FileHeader.ChannelNames = char(fh(62:61+10*FileHeader.NChannels)');
                    FileHeader.ChannelNames = reshape(FileHeader.ChannelNames,10,FileHeader.NChannels)';

                    FileHeader

                    data=int16(zeros((FileHeader.lastBlock - FileHeader.firstBlock) * 256,1));                % N of data will be changed later 

                    % read data
                    an=1;
                    for block = FileHeader.firstBlock+1 : FileHeader.lastBlock          % first data block : last data block                          
                        [data(an:an+255),count,msg] = fread(s, 256, 'int16');
                        an = an + 256;
                        fwrite(s,[hex2dec('AA'),220],'uint8')                           % confirm
                        blocks_read = blocks_read+1;
                        waitbar(double(blocks_read)/double(blocks_to_read),h)

                    end
                    data = data(1:FileHeader.NSamples);
                    data=reshape(data,FileHeader.NChannels,FileHeader.NSamples/FileHeader.NChannels)';


                    % store data on file
                    FileHeader.ScaleFactors = ones(FileHeader.NChannels,1);
                    FileHeader.Offsets = zeros(FileHeader.NChannels,1);
                    FileName = [char(FileHeader.Fileprefix),sprintf('%04d',FileHeader.FileNumber),'.mat'];
                    FileName = strrep(FileName,' ','_');
                    path='data/raw data/';
                    new_filename = inputdlg('Enter the new name for the file','File renaming');
                    FileHeader.FileName = new_filename{1};
                    save([path,new_filename{1}],'FileHeader','data')
                    fprintf('Ready : %s\n',FileName)
               end
               close(h)
            end
       case 190                                                         % setup lesen & ändern & senden
            [mat_mess,count,msg] = fread(s, len_mess, 'uint16');              % Mat mess einmal lesen
            if isempty(msg) && (mat_mess(1)== hex2dec('AAAA'))
            blocksize=mat_mess(5);
                
            [su,count,msg] = fread(s, blocksize, 'uint8');            % read setup
                    
            setup(1).FileNumber = typecast(uint8(su(1:2)),'uint16');
            setup(1).Fileprefix = char(su(3:12)');
            setup(1).TrigInt = typecast(uint8(su(13:14)),'uint16');
            setup(1).TrigDur = typecast(uint8(su(15:16)),'uint16');
            setup(1).SampFreq = typecast(uint8(su(17:18)),'uint16'); 
            setup(1).recMode = su(19);
            setup(1).Mag_set = su(20:23)';
            setup(1).ACC_set = su(24:27)';
            setup(1).IMU_set = su(28:31)';
            setup(1).ChannelNames = char(su(32:31+230)');


            disp(setup(1))

             task2=input('Change Setup (y/n) : ','s');
             if task2 == 'y'
                 setup(2)=setup(1);                    % setup 2 ist neu


                 txt2={' ';...
                     'Change what :';...
                     '   1    Fileprefix';...
                     '   2    Trigger intevall';...
                     '   3    Trigger duration';...
                     '   4    Sampling frequency';...
                     '   5    Settings for magnet';...
                     '   6    Settings for accelerometer';...
                     '   7    Settings for gyros';...                     
                     '   9    ready -> send to GaitWatch';...
                     '   10   Exit withouth updating settings';...
                    };
                task3=0;
                while task3 < 10 
                disp(char(txt2))
                 task3=input('Bitte eingeben: ');
                     switch task3
                     case 1
                         setup(2).Fileprefix='          ';
                         txt3=input('Fileprefix (max 10 characters) : ','s');
                         setup(2).Fileprefix(1:length(txt3)) = txt3;
                         setup(2).FileNumber=1;                                      % filenumber reset !!
                         
                     case 2
                         setup(2).TrigInt=input(sprintf('TrigInt in units of 1.9531 ms (%d) : ',setup(1).TrigInt));
                         
                     case 3
                         setup(2).TrigDur=input(sprintf('TrigDur in units of 1.9531 ms (%d) : ',setup(1).TrigDur));
                         
                     case 4
                         setup(2).SampFreq=input(sprintf('SampFreq (50/100/200/400) (%d) : ',setup(1).SampFreq));
                         
                     case 5    
                         for n=1:4
                             setup(2).Mag_set(n) = input(sprintf(' MagSet value %d (%d) :',n,setup(1).Mag_set(n)));
                         end
                         
                     case 6    
                         for n=1:4
                             setup(2).ACC_set(n) = input(sprintf(' ACCSet value %d (%d) :',n,setup(1).ACC_set(n)));
                         end

                     case 7    
                         for n=1:4
                             setup(2).IMU_set(n) = input(sprintf(' MagSet value %d (%d) :',n,setup(1).IMU_set(n)));
                         end

                     case 9
                            fwrite(s,[hex2dec('AA'),191],'uint8')                           % send command 191: here comes setup 
                            fwrite(s,setup(2).FileNumber,'uint16');
                            fwrite(s,setup(2).Fileprefix,'char')
                            fwrite(s,setup(2).TrigInt','uint16')
                            fwrite(s,setup(2).TrigDur','uint16')
                            fwrite(s,setup(2).SampFreq,'uint16')                            
                            fwrite(s,setup(2).recMode,'uint8')
                            fwrite(s,setup(2).Mag_set,'uint8')
                            fwrite(s,setup(2).ACC_set,'uint8')
                            fwrite(s,setup(2).IMU_set,'uint8')
                            fwrite(s,setup(2).ChannelNames,'char')
                            task3 = 10;                                   % get out
                     end
               end
             end
           end
        case 99
        OK=0;
     end
end



