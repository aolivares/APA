function [ mean_signals ] = aligned_signals( input_s1,input_s2 )
%ALIGNED_SIGNALS align two signals in time using cross-correlation and 
% carry out the mean of both.
% - Input:
%    |_ 's1,s2': input signals.
   
% - Output:
%    |_ ' mean_signals ': output signal with the same dimention than shortest.
%    It's the mean of both signals.
% -------------------------------------------------------------------------
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Veronica Torres (2): 
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
% * Last modification: 29/05/2015
% -------------------------------------------------------------------------

%   % Check the length.
%     if length(input_s1)>length(input_s2) 
        s1 = input_s1; s2 = input_s2;
%     else s1 = input_s2; s2 = input_s1;
%     end  

  % Cross correlation between both sycles.
  [acor,lag] = xcorr(s1,s2);
  
  % Selection the index where the correlation has a maximum.
  [~,I] = max(abs(acor));
  lagDiff = abs(lag(I));
  
  % Align the signals
  s1al = s1(lagDiff:end-1);
  
  % Check whether the dimensions are agree. If it's necessary, it will be
  % added samples.
  if(length(s1al)<length(s2))
      s1al=[ s1al s2(length(s1al):length(s2)-1)];
  end
  % Mean of both signals.
  mean_signals = (s1al + s2)./2;

end

