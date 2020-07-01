function [zMeas, wMeas, aMeas, zTrue, R] = measSensorReading(k, gpsLLARef)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html
global stdGyro stdAcc

% It is assumed that Euler angles are available on board for correction. Values
% avaialbe in dataset. In practicle scenario, they are calculated with 
% the help of accelerometer and magnetometer using earth gravity and earth 
% magnetic field as a reference.

stdGpsPos = [0.1, 0.1, 0.3]; %Stand. Dev. of GPS Pos Measurement, meters 
stdGpsVel = [0.05, 0.05, 0.1]; %Std. dev. of GPS Vel Measurement, meters/sec
stdEuler = [0.01, 0.01, 0.01]; %Std. dev. of euler ang. meas. error, rad/sec
R = diag([stdGpsPos, stdGpsVel, stdEuler]); %Sensor meas. noise cov. matrix


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
pos_accuracyIdx = 24;%:  velocity accuracy (north/east in m)
vel_accuracyIdx = 25;%:  velocity accuracy (north/east in m/s)
navstatIdx = 26;%:       navigation status (see navstat_to_string)
numsatsIdx = 27;%:       number of satellites tracked by primary GPS receiver
posmodeIdx = 28;%:       position mode of primary GPS receiver (see gps_mode_to_string)
velmodeIdx = 29;%:       velocity mode of primary GPS receiver (see gps_mode_to_string)
orimodeIdx = 30;%:       orientation mode of primary GPS receiver (see gps_mode_to_string)

% Load data file
fileName = strcat('000000', num2str(k, '%0.4d'), '.txt');
data = load(fileName);


gpsLLA = [data(latIdx), data(lonIdx), data(altIdx)]; %GPS reading
gpsECEF = wgslla2xyz(gpsLLA(1), gpsLLA(2), gpsLLA(3));
gpsENU = wgsxyz2enu(gpsECEF, gpsLLARef(1), gpsLLARef(2), gpsLLARef(3));
gpsNED = [gpsENU(2), gpsENU(1), -gpsENU(3)]'; 
rTrue = gpsNED;
vTrue = [data(vnIdx), data(veIdx), -data(vuIdx)]';
eulerTrue = [data(rollIdx), data(pitchIdx), data(yawIdx)]';

zTrue = [rTrue; vTrue; eulerTrue];
zMeas = zTrue + diag(R) .* randn(size(zTrue));
aTrue = [data(axIdx), data(ayIdx), data(azIdx)]';
wTrue = [data(wxIdx), data(wyIdx), data(wzIdx)]';
wMeas = wTrue + stdGyro * rand(3,1);
aMeas = aTrue + stdAcc * rand(3,1);
end