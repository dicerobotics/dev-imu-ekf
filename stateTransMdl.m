function [xPred] = stateTransMdl(xPrev, dt, aMeasPrev, wMeasPrev, nRb)%, stdAcc, stdGyro, stdDriftGyro, stdDriftAcc)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html

% State Transition Model
% x = [r, v, q, wb, ab] State Vector, (16 x 1)
% r = [rx, ry, rz]      Position, NED
% v = [vx, vy, vz]      Velocity, NED
% q = [q0, q1, q2, q3]  Body Attitude, q0 is scalar
% wb = [wbx, wby, wbz]  Angular Rate Bias, Body Frame
% ab = [abx, aby, abz]  Accelerometer Bias, Body Frame

% aMeasPrev -> Measured accelerameter reading in body frame
% wMeasPrev -> Measured gyroscope reading in body frame

global stdGyro stdAcc stdDriftDotGyro stdDriftDotAcc
global gyroOffSet accOffSet
%  gyro Offsets and accelerometer Offsets are assumed to be zero therefore
%  gyro bias and gyro drift are equal. similarly accelerometer bias and
%  accelerometer drift are also equal.
persistent wNoise
if isempty(wNoise)
    wNoise = zeros(3,1);
end
rPrev = xPrev(1:3);
vPrev = xPrev(4:6);
qPrev = xPrev(7:10);
wBiasPrev = xPrev(11:13);
aBiasPrev = xPrev(14:16);

wNoise = wNoise + stdGyro * randn(size(wMeasPrev)); %Bug Suspicion, Gyro's angular random walk (ARW)
aNoise = 0 + stdAcc * randn(size(aMeasPrev)); %Accelerometer gaussian random noise
nRb = quat2rot(qPrev); %Please test it, bug suspected

%% State Transition Vection
% Position Prediction
rPred = rPrev + vPrev * dt;
% Velocity Prediction
aGTrue = [0 0 9.80665]';    %Gravitational Acceleration, NED Frame
aGSensor = -aGTrue;         %Reason: See Important Note Below
vPred = vPrev + (nRb * (aMeasPrev - aBiasPrev) - aGSensor) * dt;
% Quaternion Prediction
qPred = (eye(4) + (dt/2) * (s2b(wMeasPrev) - s2b(wBiasPrev)))*qPrev;
% Angular Bias Prediction
wDriftPrev = wBiasPrev - gyroOffSet;
wDriftPred = eye(3) * wDriftPrev;
wBiasPred = gyroOffSet + wDriftPred;
% Angular Acceleration Prediction
aDriftPrev = aBiasPrev - accOffSet;
aDriftPred = eye(3) * aDriftPrev;
aBiasPred = accOffSet + aDriftPred;
% State Vector Prediction without process noise
xPredNoiseFree = [rPred; vPred; qPred; wBiasPred; aBiasPred];

%% Process Noise Vector
% Position Elements
wR = - nRb * aNoise * dt^2;
% Velocity Elements
wV = - nRb * aNoise * dt;
% Quaternion Elements
q0 = qPrev(1); q1 = qPrev(2); q2 = qPrev(3); q3 = qPrev(4);
qCrossMat = crossProdMat([q1, q2, q3]');
iT = [-[q1, q2, q3]; ...
        q0*eye(3)+qCrossMat];
wQ = -(dt/2) * iT * wNoise;
% Angular Bias(Drift) Elements
wWDriftDot = stdDriftDotGyro * randn(3,1) * dt;
% Angular Acceleration Elements
wADriftDot = stdDriftDotAcc * randn(3,1) * dt;

% Process Noise Vector
wProcessNoiseVec = [wR; wV; wQ; wWDriftDot; wADriftDot];

%%
% State vector prediction with process noise.
xPred = xPredNoiseFree;% + wProcessNoiseVec;

%%
% Important Note: We need to multiply Gravity Vector aG by negative one.
% Along the axis pointed in the direction of gravity, the accelerometer measures -1 [g]. 
% This is due to the proof-mass being pulled in the direction of gravity, which is equivalent 
% to a deceleration of 1 [g] in the absence of gravity.
% Ref: https://openimu.readthedocs.io/en/latest/algorithms/MeasurementVector.html
end

