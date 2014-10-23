% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% |||||||||||||||||||||| GAIT ANALYSIS ROUTINE ||||||||||||||||||||||||||||
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% -------------------------------------------------------------------------
% 1) Load, extract, calibrate data and compute orientation.
% -------------------------------------------------------------------------
main;

% -------------------------------------------------------------------------
% 2) Compute gait parameters.
% -------------------------------------------------------------------------
% Check if all the legs' segments were selected.
aux1 = strfind(S(Selection),'right shank');
aux1 = aux1(~cellfun('isempty',aux1));
aux2 = strfind(S(Selection),'left shank');
aux2 = aux2(~cellfun('isempty',aux2));
aux3 = strfind(S(Selection),'right thigh');
aux3 = aux3(~cellfun('isempty',aux3));
aux4 = strfind(S(Selection),'left thigh');
aux4 = aux4(~cellfun('isempty',aux4));

% If so, the gait analysis routine is called.
if ~isempty(aux1) && ~isempty(aux2) && ~isempty(aux3) && ~isempty(aux4)
% HERE SHOULD GO THE CALL TO THE FUNCTION WHICH COMPUTES THE GAIT
% PARAMETERS. THE FOLLOWING SIGNALS ARE AVAILABLE FOR SUCH COMPUTATION:

% * Acceleration & Angular rate:
% - 'a_X_right_shank_1_C' -> Calibrated acceleration (X axis, right shank).
% - 'a_Z_right_shank_1_C' -> Calibrated acceleration (Z axis, right shank).
% - 'g_Y_right_shank_1_C' -> Calibrated angular rate (Y axis, right shank).
% - 'a_X_left_shank_1_C' -> Calibrated acceleration (X axis, left shank).
% - 'a_Z_left_shank_1_C' -> Calibrated acceleration (Z axis, left shank).
% - 'g_Y_left_shank_1_C' -> Calibrated angular rate (Y axis, left shank).
% - 'a_X_right_thigh_1_C' -> Calibrated acceleration (X axis, right thigh).
% - 'a_Z_right_thigh_1_C' -> Calibrated acceleration (Z axis, right thigh).
% - 'g_Y_right_thigh_1_C' -> Calibrated angular rate (Y axis, right thigh).
% - 'a_X_left_thigh_1_C' -> Calibrated acceleration (X axis, left thigh).
% - 'a_Z_left_thigh_1_C' -> Calibrated acceleration (Z axis, left thigh).
% - 'g_Y_left_thigh_1_C' -> Calibrated angular rate (Y axis, left thigh).

% * Orientation angles:
% - 'pitch_KF_right_shank' -> pitch of right shank.
% - 'pitch_KF_left_shank'  -> pitch of left shank.
% - 'pitch_KF_right_thigh' -> pitch of right thigh.
% - 'pitch_KF_left_thigh'  -> pitch of left thigh.
    
else
    msgbox('Error: All segments of the leg must be selected.')
end