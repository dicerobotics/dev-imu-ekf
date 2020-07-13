function [H] = symObservationMdl(xPred)
% Script Writer:	Rana Zain Idrees
% Association:      UET, Lahore
% Date Created:     July 7th, 2020
% email:            rzi_93@hotmail.com
% Advisor:          Dr. K. M. Hassan
% Code Ref:         MATLAB Document
% Modified by:      Awais Arshad
% Date Modified:	July 7th, 2020

obsMdlR = eye(3);
obsMdlV = eye(3);

% syms q0 q1 q2 q3 r p y 
% r = atan2( 2 * ( q2 * q3 + q0 * q1 ), q0^2 - q1^2 - q2^2 + q3^2); % roll 
% p = - asin( 2 * ( q1 * q3 - q0 * q2 )); % pitch 
% y = atan2( 2 * ( q1 * q2 + q0 * q3 ), q0^2 + q1^2 - q2^2 - q3^2); % yaw 
% q = [q0; q1; q2; q3]; % 4x1 
% eulerAngles = [r; p; y]; % 3x1
% symHQuat = jacobian(eulerAngles,q); %Hardcoded below for comput. effi.
qPred = xPred(7:10); q0 = qPred(1);  q1 = qPred(2); q2 = qPred(3); q3 = qPred(4);
% HQuat = subs(symHQuat); Primary solution
% obsMdlEuler = double(HQuat);


% Modified solution - Hardcoded version of syms method
HQuat(1,1) = (((2*imag(q0) + 2*real(q1))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i)) + ((2*imag(q1) - 2*real(q0))*(imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i)))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)*(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)/((real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2 + (imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i))^2);
HQuat(1,2) =  -(((2*imag(q1) - 2*real(q0))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i)) - ((2*imag(q0) + 2*real(q1))*(imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i)))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)*(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)/((real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2 + (imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i))^2);
HQuat(1,3) = -(((2*imag(q2) - 2*real(q3))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i)) - ((2*imag(q3) + 2*real(q2))*(imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i)))/(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)*(real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2)/((real(q0^2) - real(q1^2) - real(q2^2) + real(q3^2) + real(q0*q1*2i) + real(q2*q3*2i))^2 + (imag(q0^2) - imag(q1^2) - imag(q2^2) + imag(q3^2) + imag(q0*q1*2i) + imag(q2*q3*2i))^2);
HQuat(1,4) = -(((2*imag(q3) - 2*real(q0))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) - ((2*imag(q0) + 2*real(q3))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);
HQuat(2,1) = (2*q2)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
HQuat(2,2) = -(2*q3)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
HQuat(2,3) = -(((2*imag(q3) - 2*real(q0))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) - ((2*imag(q0) + 2*real(q3))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);
HQuat(2,4) = -(2*q1)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
HQuat(3,1) = (((2*imag(q0) + 2*real(q3))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) + ((2*imag(q3) - 2*real(q0))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);
HQuat(3,2) = (((2*imag(q1) + 2*real(q2))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) + ((2*imag(q2) - 2*real(q1))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);
HQuat(3,3) =  -(((2*imag(q3) - 2*real(q0))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) - ((2*imag(q0) + 2*real(q3))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);
HQuat(3,4) = -(((2*imag(q3) - 2*real(q0))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i)) - ((2*imag(q0) + 2*real(q3))*(imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i)))/(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)*(real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2)/((real(q0^2) + real(q1^2) - real(q2^2) - real(q3^2) + real(q0*q3*2i) + real(q1*q2*2i))^2 + (imag(q0^2) + imag(q1^2) - imag(q2^2) - imag(q3^2) + imag(q0*q3*2i) + imag(q1*q2*2i))^2);

obsMdlEuler = HQuat;
H = [obsMdlR, zeros(3,13); ... 
    zeros(3,3), obsMdlV, zeros(3,10); ... 
    zeros(3,6), obsMdlEuler, zeros(3,6)];

end

