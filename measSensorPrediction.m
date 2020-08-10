function h = measSensorPrediction(xPred)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html

% xPred = [rPred; vPred; qPred; wBiasPred; aBiasPred];
rPred = xPred(1:3);
vPred = xPred(4:6);
qPred = xPred(7:10);
% eulerPred = quat2euler(qPred);
[yaw, pitch, roll] = quat2angle(qPred', 'ZYX');
eulerBody2NED = [roll, pitch, yaw]';
eulerPred = eulerBody2NED;
h = [rPred; vPred; eulerPred];
end
