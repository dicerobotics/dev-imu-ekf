function [zMeas, wMeas, aMeas, zTrue, R] = measSensorReading(k, gpsLLARef)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         See function main();

global stdGyro stdAcc stdDriftDotGyro stdDriftDotAcc
global stdGpsPos stdGpsVel stdEuler
global gyroOffSet accOffSet
persistent wDriftTruePre aDriftTruePre
if (isempty(wDriftTruePre) && isempty(aDriftTruePre))
    wDriftTruePre = zeros(3,1);
    aDriftTruePre = zeros(3,1);
end

% It is assumed that Euler angles are available on board for correction. Values
% avaialbe in dataset. In practical scenario, they are calculated with 
% the help of accelerometer and magnetometer using earth gravity and earth 
% magnetic field as a reference.


%Define Indices accroding to data files
latIdx = 1;%:   latitude of the oxts-unit (deg)
lonIdx = 2;%:   longitude of the oxts-unit (deg)
altIdx = 3;%:   altitude of the oxts-unit (m)

rollIdx = 4;%:  roll angle (rad),    0 = level, positive = left side up,      range: -pi   .. +pi
pitchIdx = 5;%: pitch angle (rad),   0 = level, positive = front down,        range: -pi/2 .. +pi/2
yawIdx = 6;%:   heading (rad),       0 = east,  positive = counter clockwise, range: -pi   .. +pi

vnIdx = 7;%:    velocity towards north (m/s)
veIdx = 8;%:    velocity towards east (m/s)
vfIdx = 9;%:    forward velocity, i.e. parallel to earth-surface (m/s)
vlIdx = 10;%:    leftward velocity, i.e. parallel to earth-surface (m/s)
vuIdx = 11;%:    upward velocity, i.e. perpendicular to earth-surface (m/s)

axIdx = 12;%:    acceleration in x, i.e. in direction of vehicle front (m/s^2)
ayIdx = 13;%:    acceleration in y, i.e. in direction of vehicle left (m/s^2)
azIdx = 14;%:    acceleration in z, i.e. in direction of vehicle top (m/s^2)

afIdx = 15;%:    forward acceleration (m/s^2)
alIdx = 16;%:    leftward acceleration (m/s^2)
auIdx = 17;%:    upward acceleration (m/s^2)

wxIdx = 18;%:    angular rate around x (rad/s)
wyIdx = 19;%:    angular rate around y (rad/s)
wzIdx = 20;%:    angular rate around z (rad/s)

wfIdx = 21;%:    angular rate around forward axis (rad/s)
wlIdx = 22;%:    angular rate around leftward axis (rad/s)
wuIdx = 23;%:    angular rate around upward axis (rad/s)

posAccuracyIdx = 24;%:  velocity accuracy (north/east in m)
velAccuracyIdx = 25;%:  velocity accuracy (north/east in m/s)

navstatIdx = 26;%:       navigation status (see navstat_to_string)
numsatsIdx = 27;%:       number of satellites tracked by primary GPS receiver
posmodeIdx = 28;%:       position mode of primary GPS receiver (see gps_mode_to_string)
velmodeIdx = 29;%:       velocity mode of primary GPS receiver (see gps_mode_to_string)
orimodeIdx = 30;%:       orientation mode of primary GPS receiver (see gps_mode_to_string)

%%
% Load data file
fileName = strcat('000000', num2str(k, '%0.4d'), '.txt');
data = load(fileName);

% ENU, NED, Sensor, and Body Frame Relationships
% ENU: x->East, y-> North, z->Up (Zenith)
% NED:  x->North, y-> East, z->Down (Nadir)
% Sensor:  x->Forward, y-> Left, z->Top
% Body:   x->Forward, y->Right, z->Down
eulerENU = [data(rollIdx), data(pitchIdx), data(yawIdx)]'; %Sensor w.r.t. ENU
sRb = euler2rot([pi, 0, 0]');       %Rotation matrix of Body frame represented in Sensor frame
eRs = euler2rot(eulerENU);          %Rotation matrix of ENU frame represented in Body frame
% eRs = sRe';                         %Rotation matrix for Sensor frame represented in ENU frame
nRe = [0, 1, 0; 1, 0, 0; 0, 0, -1]; %Rotation matrix for ENU frame represented in NED frame
nRs = nRe * eRs;                    %Rotation matrix for Sensor frame represented in NED frame
nRb = nRe * eRs * sRb;              %Rotation matrix for Body frame represented in NED frame
% disp('nRb New inside measSensorReading'); disp(nRb);


% Measurement Noise Covariance Matrix
R = diag([stdGpsPos, stdGpsVel, stdEuler].^2); %Sensor meas. noise cov. matrix
% R = diag([0.5 0.5 0.5 0.0139 0.0139 0.0139 0 0 0]);
% R = zeros(9,9);
% R = diag([data(posAccuracyIdx)*ones(1,3), data(velAccuracyIdx)*zeros(1,3), stdEuler].^2); %Sensor meas. noise cov. matrix
%%
% Position Readout
gpsPosLLA = [data(latIdx), data(lonIdx), data(altIdx)]; %GPS reading
gpsPosECEF = wgslla2xyz(gpsPosLLA(1), gpsPosLLA(2), gpsPosLLA(3));
gpsPosENU = wgsxyz2enu(gpsPosECEF, gpsLLARef(1), gpsLLARef(2), gpsLLARef(3));
gpsPosNED = [gpsPosENU(2), gpsPosENU(1), -gpsPosENU(3)]'; 
rTrue = gpsPosNED;

% Velocity Readout
vNED = [data(vnIdx), data(veIdx), -data(vuIdx)]'; %NED Frame, Assumption: upword vector from earth surface is parallel to wgs84 upword vector at referece LLA
vTrue = vNED;

% Euler Readout
eulerNED = rot2euler(nRb);	%WTF %Body w.r.t. NED, Rotation Sequence: ZYX
eulerTrue = eulerNED;   %Consider constraining euler angles here

zTrue = [rTrue; vTrue; eulerTrue];
zErr = sqrt(diag(R)) .* randn(size(zTrue));
zMeas = zTrue + zErr;

aSensKitti = [data(axIdx), data(ayIdx), data(azIdx)]'; %SensKitti: x->Forward, y->Left, z->Up
% aSensKitti(3) = -aSensKitti(3); %Reason: See Important Note Below
aBody = [aSensKitti(1), -aSensKitti(2), -aSensKitti(3)]'; %Body: x->Forward, y->Right, z->Down
aTrue = aBody;

wSensKitti = [data(wxIdx), data(wyIdx), data(wzIdx)]'; %SensKitti: x->Forward, y->Left, z->Up
wBody = [wSensKitti(1), -wSensKitti(2), -wSensKitti(3)]'; %Body: x->Forward, y->Right, z->Down
wTrue = wBody;

wDriftTrue = wDriftTruePre + stdDriftDotGyro * randn(3,1);
aDriftTrue = aDriftTruePre + stdDriftDotAcc * randn(3,1);
wThermalNoiseTrue = stdGyro * randn(3,1);
aThermalNoiseTrue = stdAcc * randn(3,1);
wMeas = wTrue + gyroOffSet + wDriftTrue +  wThermalNoiseTrue;
aMeas = aTrue + accOffSet + aDriftTrue + aThermalNoiseTrue;
wDriftTruePre = wDriftTrue; aDriftTruePre = aDriftTrue;
end