function [coordenates] = getCoordenates( signal_x,signal_y, fig_title )

% FUNCTION GETCOORDENATES extracts the coordenates of the data points shown  
% in a figure and selected using the data cursor tool.
%
% Input parameters:
% |_ 'signal': The signal to be plotted.
% |_ 'fig_title': Title which is shown in the figure. 
%
% Output parameters:
% |_ 'coordenates': The selected coordenates.
%
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
% * Last modification: 05/06/2015
% -------------------------------------------------------------------------
% Create figure displaying the signal.
    fig = figure; 
    plot(signal_x,signal_y,'+')
    set(gcf,'position',[40 80 1600 900])
    title(fig_title);

% The data cursor mode is activated. The execution of the routine does not 
    % continue until the user clicks the "Continue execution" button. 
    set(gcf,'position',[40 80 1600 900])
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')
    uicontrol('Style', 'pushbutton','String', 'Continue execution',...
        'Position', [720 20 150 20],'Callback', 'uiresume(gcbf)');
    uiwait(fig);

    % With 'Position' option we can obtain the coordenates of the point
    % selected.
    c_info = getCursorInfo(dcm_obj);
    indexes = zeros(1,length(c_info));
    for i = 1:length(c_info)
        coordenates (i,:) = c_info(i).Position;
    end

end

