function [ peaks_index ] = findMaxPeaks( signal, initcross, finalcross, sing)
%FINDMAXPEAKS find the maximum peak (negative or positive) in a 
%interval specified. If any peak is found, it will be calculated the
%maximum orminimum value in each case.
% - Input:
%    |_ 'signal': input signal.
%    |_ 'initcross': beginning of the   the peak will be looked
%    for.
%    |_ 'finalcross': final of the interval where the peak will be looked
%    for.
%    |_ 'sing': if sing is==1 (we will look for a positive peak) and if it
%    is ==2 (we will look for a negative peak).
%    
% - Output:
%    |_ ' peaks_index ': index where the peak have been found.
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
% * Last modification: 22/05/2015
% -------------------------------------------------------------------------

% Differenciate when we serch for a positive or negative peak.

    if(sing==1)% Search for a positive peak

         % Find all peaks in each interval.
        [pos_peak_values, pos_peak_locations] = findpeaks(...
                                                signal(initcross:finalcross));

         % Check if there are positive peaks detected.
         if isempty(pos_peak_values)
                [pos_peak_values, pos_peak_locations] = max( signal(...
                                                    initcross:finalcross)); 
                peaks_index = pos_peak_locations + initcross - 1;
          else

         % Store the index of the longest positive peak.                                      
                peaks_index = find(signal(initcross:finalcross)== max(...
                                          pos_peak_values), 1) + initcross - 1;
         end  
         
    elseif(sing==2)% Search for a negative peak.
        
        % Find all peaks in each interval.
        [neg_peak_values, neg_peak_locations] = findpeaks(...
                                                -signal(initcross:finalcross));

         % Check if there are negative peaks detected.
         if isempty(neg_peak_values)
           [neg_peak_values, neg_peak_locations] = min(...
                                            signal(initcross:finalcross));
           peaks_index = neg_peak_locations + initcross -1;
                                        
         else
            % Store the index of the longest positive peak.                                      
            peaks_index = find(signal(initcross:finalcross)== -max(...
                                      neg_peak_values), 1) + initcross - 1;                                    
         end


     end
end

