function [argOut] = repairQuaternion(argIn)
% Script Writer:	Rana Zain Idrees
% Association:      UET, Lahore
% Date Created:     July 7th, 2020
% email:            rzi_93@hotmail.com
% Advisor:          Dr. K. M. Hassan
% Code Ref:         Unknown
% Modified by:      Awais Arshad
% Date Modified:	July 8th, 2020

%% 
% Quat Extraction
% Quat = [scalar; vec3x1];
if (length(argIn)>4)
    quatIdx = 7:10; %considering whole state vector as an input
else
    quatIdx = 1:4;  %considering raw quaternion as an input
end

qparts = argIn(quatIdx);

% Quat Normallization 
n = sqrt(sum(qparts.^2));
qparts = qparts./n;

% Correction
argOut = argIn; %To cover the case of whole state vector as an input
if qparts(1) < 0
    argOut(quatIdx) = -qparts;
else
    argOut(quatIdx) = qparts;
end

end

