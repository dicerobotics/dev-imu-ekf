function [H] = obsMdl(xPred)
% Script Writer:	Rana Zain Idrees
% Association:      UET, Lahore
% Date Created:     July 7th, 2020
% email:            rzi_93@hotmail.com
% Advisor:          Dr. Khalid Mehmood Hassan
% Code Ref:         MATLAB Documents
% Modified by:      Awais Arshad
% Date Modified:	Aug 5th, 2020

obsMdlR = eye(3);
obsMdlV = eye(3);

persistent hndlHQuat
if isempty(hndlHQuat)
    syms q0 q1 q2 q3 r p y 
    r = atan2( 2 * ( q2 * q3 + q0 * q1 ), q0^2 - q1^2 - q2^2 + q3^2); % roll 
    p = - asin( 2 * ( q1 * q3 - q0 * q2 )); % pitch 
    y = atan2( 2 * ( q1 * q2 + q0 * q3 ), q0^2 + q1^2 - q2^2 - q3^2); % yaw 
    q = [q0;  q1;  q2;  q3]; % 4x1 
    eulerAngles = [r; p; y]; % 3x1
    symHQuat = jacobian(eulerAngles,q);
    hndlHQuat = matlabFunction(symHQuat);
end

qPred = xPred(7:10); q0 = qPred(1);  q1 = qPred(2); q2 = qPred(3); q3 = qPred(4);
HQuat = hndlHQuat(q0, q1, q2, q3);
obsMdlEuler = HQuat;

H = [obsMdlR, zeros(3,13); ... 
    zeros(3,3), obsMdlV, zeros(3,10); ... 
    zeros(3,6), obsMdlEuler, zeros(3,6)];

end

