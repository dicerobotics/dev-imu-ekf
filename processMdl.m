function [F, Q] = processMdl(xPrev, dt, aMeasPrev, wMeasPrev)%, stdAcc, stdGyro, stdDriftGyro, stdDriftAcc)

% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html

global stdGyro stdAcc stdDriftDotGyro stdDriftDotAcc
qPrev = xPrev(7:10);
wBiasPrev = xPrev(11:13);
aBiasPrev = xPrev(14:16);
aPrev = aMeasPrev - aBiasPrev;

% aG = [0 0 9.8]'; % Gravitational Acceleration, NED Frame
% nRb = quat2rot(qPrev)'; bRn = nRb';
% aPrev = bRn*(nRb * (aMeasPrev - aBiasPrev) - aG);
q0 = qPrev(1); q1 = qPrev(2); q2 = qPrev(3); q3 = qPrev(4);
Qf = [q1, q0, -q3, q2; ... 
      q2, q3, q0, -q1; ...
      q3, -q2, q1, q0];

aCrossMat = crossProdMat(aPrev);
pvpq = 2 * Qf * [0, aPrev'; ... 
                aPrev, -aCrossMat];

qCrossMat = [0 -q3 q2; ...
             q3, 0 -q1; ...
            -q2, q1, 0];

qV = [qPrev(2), qPrev(3), qPrev(4)];
iT = [-qV; ...
       qPrev(1)*eye(3)+qCrossMat];

bigOmg = (s2b(wMeasPrev) - s2b(wBiasPrev));
nRb = quat2rot(qPrev)';
z3 = zeros(3,3); i3 = eye(3); z34 = zeros(3,4); z43 = zeros(4,3); 
            
f = [z3, i3, z34, z3, z3; ...
    z3, z3, pvpq, z3, -nRb; ...
    z43, z43, (1/2)*bigOmg, -(1/2)*iT, z43; ...
    z3, z3, z34, z3, z3; ...
    z3, z3, z34, z3, z3];
F = eye(16) + f * dt;

covR = (stdAcc * dt^2)^2 * eye(3);
covV = (stdAcc * dt)^2 * eye(3);
covQ = (stdGyro*dt/2)^2 * [1-q0^2,-q0*q1,-q0*q2,-q0*q3; ... 
                           -q0*q1,1-q1^2,-q1*q2,-q1*q3; ...
                           -q0*q2,-q1*q2,1-q2^2,-q2*q3; ... 
                           -q0*q3,-q1*q3,-q2*q3,1-q3^2];
covWBias = (stdDriftDotGyro * dt )^2 * eye(3);
covABias = (stdDriftDotAcc * dt )^2 * eye(3);

Q = [covR, z3, z34, z3, z3; ... 
    z3, covV, z34, z3, z3; ... 
    z43, z43, covQ, z43, z43; ... 
    z3, z3, z34, covWBias, z3; ... 
    z3, z3, z34, z3, covABias]; 

%Q -> Process noise covariance matrix
%P -> Process estimation error covariance matrix
%R -> Measurement Noise Covariance Matrix
end
