function H = observationMdl(xPred)

% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         Unknown

qPred = xPred(7:10);
obsMdlR = eye(3);
obsMdlV = eye(3);

q0 = qPred(1); q1 = qPred(2); q2 = qPred(3); q3 = qPred(4);

u = 2*(q2*q3+q0*q1); 
pupq0 = 2*q1; pupq1 = 2*q0; pupq2 = 2*q3; pupq3 = 2*q2;

v = q0^2-q1^2-q2^2+q3^2;
pvpq0 = 2*q0; pvpq1 = -2*q1; pvpq2 = -2*q2; pvpq3 = 2*q3;

w = 2*(q1*q3-q0*q2);
pwpq0 = -2*q2; pwpq1 = 2*q3; pwpq2 = -2*q0; pwpq3 = 2*q1;

x = 2*(q1*q2+q0*q3); 
pxpq0 = 2 * q3; pxpq1 = 2 * q2; pxpq2 = 2 * q1; pxpq3 = 2 * q0;

y = (q0^2+q1^2-q2^2-q3^2);
pypq0 = 2 * (q0); pypq1 = 2 * (q1); pypq2 = 2 * (-q2); pypq3 = 2 * (-q3);

obsMdlRoll = (v^2/(v^2+u^2)) * [v*pupq0-u*pvpq0, v*pupq1-u*pvpq1, v*pupq2-u*pvpq2, v*pupq3-u*pvpq3];
obsMdlPitch = (-1/sqrt(1-w^2)) * [pwpq0, pwpq1, pwpq2, pwpq3];
obsMdlYaw = (y^2/(y^2+x^2)) * [y*pxpq0-x*pypq0, y*pxpq1-x*pypq1, y*pxpq2-x*pypq2, y*pxpq3-x*pypq3];

obsMdlEuler = [obsMdlRoll; obsMdlPitch; obsMdlYaw];

% x = [r, v, q, wb, ab] State Vector, (16 x 1)
% h = [r, v, euler]
H = [obsMdlR, zeros(3,13); ... 
    zeros(3,3), obsMdlV, zeros(3,10); ... 
    zeros(3,6), obsMdlEuler, zeros(3,6)];
end