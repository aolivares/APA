function F = eofGKF(p)

global gyro;
global acc;
global frec;
global true_angle;
global marker;

gw = gwLibrary;

alpha1 = p(1);
alpha2 = p(2);
beta1 = p(3);
beta2 = p(4);

angle_DKF = gw.fusionGKF(gyro,acc,frec,var(acc),var(acc),...
    var(gyro),alpha1,alpha2,beta1,beta2,acc(1),marker);

F = sqrt(mean((true_angle-angle_DKF').^2));