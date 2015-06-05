% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% --------------------------- APA STADISTICS ------------------------------
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% * Project name: Comparison of Posturographic Body-sway Measurements with 
%                 Accelerometric Data.
%
% * Authors:      - Prof. Dr. Med. Kai Boetzel (1): 
%                   |_ kai.boetzel@med.uni-muenchen.de 
%                 - Veronica  Torres (2): 
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

% -------------------------------------------------------------------------
% 0) Clear workspace.
% -------------------------------------------------------------------------
clear all; close all; clc;
filepath =  '../../data/APA Parameters/FPvsGW';

load(fullfile(filepath,'ES39'));
PCA_x_1 = peaks_APA_PCA_x;
PCA_y_1 = peaks_APA_PCA_y;

load(fullfile(filepath,'MM57'));
PCA_x_2 = peaks_APA_PCA_x;
PCA_y_2 = peaks_APA_PCA_y;

load(fullfile(filepath,'RS46'));
PCA_x_3 = peaks_APA_PCA_x;
PCA_y_3 = peaks_APA_PCA_y;

load(fullfile(filepath,'SW47'));
PCA_x_4 = peaks_APA_PCA_x;
PCA_y_4 = peaks_APA_PCA_y;

load(fullfile(filepath,'WS42'));
PCA_x_5 = peaks_APA_PCA_x;
PCA_y_5 = peaks_APA_PCA_y;

% Compare the components
PCA_x_first = [PCA_x_1(1,:); PCA_x_2(1,:);PCA_x_3(1,:);PCA_x_4(1,:);PCA_x_5(1,:)];
PCA_x_second = [PCA_x_1(2,:); PCA_x_2(2,:);PCA_x_3(2,:);PCA_x_4(2,:);PCA_x_5(2,:)];

PCA_y_first = [PCA_y_1(1,:); PCA_y_2(1,:);PCA_y_3(1,:);PCA_y_4(1,:);PCA_y_5(1,:)];
PCA_y_second = [PCA_y_1(2,:); PCA_y_2(2,:);PCA_y_3(2,:);PCA_y_4(2,:);PCA_y_5(2,:)];

% Correlation

[corr_x_first, ~] = corr(PCA_x_first(:,1),PCA_x_first(:,2));
[corr_x_second, ~] = corr(PCA_x_second(:,1),PCA_x_second(:,2));

[corr_y_first, ~] = corr(PCA_y_first(:,1),PCA_y_first(:,2));
[corr_y_second, ~] = corr(PCA_y_second(:,1),PCA_y_second(:,2));




