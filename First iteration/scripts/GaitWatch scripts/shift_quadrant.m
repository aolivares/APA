yaw_mag_kf = atan2d(hx_n_KF,hy_n_KF);
yaw_mag_kf = gw.correct_yaw_quad_shifts(yaw_mag_kf,'deg');

correct_again = 1;
while correct_again == 1
    correct_again = 0;
    for i = 2: length(yaw_mag_kf)
        if abs(yaw_mag_kf(i)-yaw_mag_kf(i-1)) > limit
            correct_again = 1;
        end
    end
    if correct_again == 1
        yaw_mag_kf = gw.correct_yaw_quad_shifts(yaw_mag_kf,'deg');
    end
end
uncorrected_angle = yaw_mag_kf;

limit = 250;
limit2 = 361;
correction = 0;
rectification =  abs(min(uncorrected_angle));
if min(uncorrected_angle) < 0
    corrected_angle = uncorrected_angle + rectification;
    correction = 1;
else
    corrected_angle = uncorrected_angle;
end

figure
plot(uncorrected_angle)
hold on
plot(corrected_angle,'r')

corrected_angle2 = corrected_angle;
for i = 2: length(corrected_angle)
    distance = abs(corrected_angle(i) - corrected_angle(i-1));
    if  distance > limit
        if corrected_angle(i) > corrected_angle(i-1) 
            corrected_angle(i) = corrected_angle(i) - 360;
        elseif corrected_angle(i) < corrected_angle(i-1) 
            corrected_angle(i) = corrected_angle(i) + 360;
        end
    end
%     elseif distance > limit2
%         if corrected_angle(i) > corrected_angle(i-1) 
%             corrected_angle(i) = corrected_angle(i) - 720;
%         elseif corrected_angle(i) < corrected_angle(i-1) 
%             corrected_angle(i) = corrected_angle(i) + 720;
%         end
%     end  
end

figure
plot(corrected_angle2)
hold on
plot(corrected_angle,'r')

if correction == 1
    corrected_angle = corrected_angle - abs(min(uncorrected_angle));
end