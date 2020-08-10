function [] = genSyntheticData()
clc; clear; close all; 
disp('Synthetic Data Generation in Progress...');
addpath('./subModules');

startTime = datetime; %current date and time
sampleTime = 0.01; %seconds
duration = 15; %seconds
sampleTimeSeries = 0:sampleTime:duration; %seconds, First sample at zero

genTimeStamps(startTime, sampleTimeSeries);
genDataFiles(sampleTimeSeries);
end


%%
function [] = genDataFiles(sampleTimeSeries)
disp('Generating data files...');
llaRef = [49.008644826538 8.3981039999565 112.99059295654]; %Must match referece LLA in main file.
dataFilesPath = './dataSet/oxtsSynthetic/data/';
if ~exist(dataFilesPath, 'dir')
    mkdir('./dataSet/oxtsSynthetic/', 'data');
end
addpath(dataFilesPath);

accFreq = [1/10, 1/15, 1/30]'; %Hz
accAmp = 1; %m/s^2
accNED_GFree = accAmp * cos(2*pi*(accFreq.*sampleTimeSeries)); %Gravity Free
velNED= (accAmp./(2*pi.*accFreq)) .* sin(2*pi*(accFreq.*sampleTimeSeries));%Integration for acceleration
posNED = -(accAmp./(2*pi.*accFreq).^2) .* cos(2*pi*(accFreq.*sampleTimeSeries));%Integration for velocity
aGNED = [0 0 9.80665]';    %Gravitational Acceleration, NED Frame

wFreq = [1/10, 1/15, 1/5]';
wAmp = 1*[0 0 pi/4]';  %Keep Roll pitch zero. For gravitation acceleration.
wNED = wAmp .* cos(2*pi*(wFreq.*sampleTimeSeries)); % Gyroscopic signal varies as a sinosoid.
angleNed2Body = (wAmp./(2*pi.*wFreq)) .* sin(2*pi*(wFreq.*sampleTimeSeries)); %Angular Displacement
% figure(); plot(sampleTimeSeries, angleNed2Body(3,:));
% Initial State
numSamples = size(angleNed2Body,2);
bRn = angle2dcm(angleNed2Body(3,:), angleNed2Body(2,:), angleNed2Body(1,:), 'ZYX');

for n=1:numSamples
    accBodyGFree = bRn(:,:,n) * accNED_GFree(:,n); %Gravity Free Acceleleration
    aGBody = bRn(:,:,n) * aGNED;
    accBody = accBodyGFree - aGBody; %nevative sign to compensation sensor measurement principle
    
    wBody = -wNED(:,n); %I am not sure about this statement myself. But it's important.
    
    nRb = bRn(:,:,n)';
    nQb = dcm2quat(nRb);

    
    xN = [posNED(:,n); velNED(:,n); nQb']; %Nth state vector
    oxtsData = makeOxtsData(xN, llaRef, accBody, wBody);
    
    fileName = strcat(num2str(n-1, '%0.10d'), '.txt'); fileIDData = fopen([dataFilesPath, fileName], 'wt' );
    formatSpec = '%015f %015f %015f %08f %08f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %015f %1.0f %1.0f %1.0f %1.0f %1.0f';
    fprintf(fileIDData, formatSpec,oxtsData); fclose(fileIDData);
end

disp('Data files generation completed successfully');
disp('Synthetic data generation task completed successfully');
end

function [dataArray] = makeOxtsData(xCurr, llaRef, aBody, wBody)
% Read Arguments
posNED = xCurr(1:3,1);
velNED = xCurr(4:6,1);
nQb = xCurr(7:10,1);

% Calculate Data Values
wgs84 = referenceEllipsoid('wgs84');
[lat,lon,alt] = ned2geodetic(posNED(1),posNED(2),posNED(3),llaRef(1),llaRef(2),llaRef(3),wgs84);
posLLA = [lat,lon,alt]';

% nRb = quat2rot(quatNEDBody); %rotation matrix for body frame represented in euler frame
nRb = quat2dcm(nQb');
nRe = [0, 1, 0; 1, 0, 0; 0, 0, -1]; %Rotation matrix for ENU frame represented in NED frame
eRn = nRe';% eRb = eRn * nRb;
sRb = angle2dcm(0, 0, pi, 'ZYX'); bRs = sRb;
eRs = eRn * nRb * bRs; sRe = eRs';

[yaw, pitch, roll] = dcm2angle(sRe, 'ZYX');
eulerENU2Sensor = fliplr([yaw, pitch, roll])';

% Make Data Array
latIdx = 1; lonIdx = 2; altIdx = 3;
rollIdx = 4; pitchIdx = 5; yawIdx = 6;
vnIdx = 7; veIdx = 8; vuIdx = 11;
axIdx = 12; ayIdx = 13; azIdx = 14;
wxIdx = 18; wyIdx = 19; wzIdx = 20;

dataArray = -1 * ones(1,30); %Fake Initialization

dataArray(1, latIdx) = posLLA(1);
dataArray(1, lonIdx) = posLLA(2);
dataArray(1, altIdx) = posLLA(3);

dataArray(1, rollIdx) = eulerENU2Sensor(1);
dataArray(1, pitchIdx) = eulerENU2Sensor(2);
dataArray(1, yawIdx) = eulerENU2Sensor(3);

dataArray(1, vnIdx) = velNED(1);
dataArray(1, veIdx) = velNED(2);
dataArray(1, vuIdx) = -velNED(3);

dataArray(1, axIdx) = aBody(1);
dataArray(1, ayIdx) = -aBody(2);
dataArray(1, azIdx) = -aBody(3);

dataArray(1, wxIdx) = wBody(1);
dataArray(1, wyIdx) = -wBody(2);
dataArray(1, wzIdx) = -wBody(3);
end




%%
function [] = genTimeStamps(startTime, sampleTimeSeries)
disp('Generating time stamps...');
fileName = 'timestamps';

timeStampFilePath = './dataSet/oxtsSynthetic/';
if ~exist(timeStampFilePath, 'dir')
    mkdir('dataSet', 'oxtsSynthetic');
end
addpath(timeStampFilePath);

fileIDTimeStamps = fopen([timeStampFilePath, fileName, '.txt'], 'wt' );

tStart = startTime;
% sampleTimeSeries = 0:sampleTime:duration;
for dt = sampleTimeSeries
  t = tStart + seconds(dt);
  fprintf(fileIDTimeStamps, string(t, 'yyyy-MM-dd HH:mm:ss.SSSSSSSSS'));
  fprintf(fileIDTimeStamps, '\n');
end
fclose(fileIDTimeStamps);
disp('Time stamps generation completed successfully');
end



function [xProp] = stateProp(xPrev, dt, aBody, wBody)
rPrev = xPrev(1:3);
vPrev = xPrev(4:6);
qPrev = xPrev(7:10);
% nRb = quat2rot(qPrev); %Please test it, bug suspected
nRb = quat2dcm(qPrev');
%% State Transition Vection

% Velocity propagation
aGNED = [0 0 9.80665]';    %Gravitational Acceleration, NED Frame
aG_NEDSens = -aGNED;         %Reason: See Important Note Below, Roll/Pitch Must be Zero
vPred = vPrev + (nRb * aBody - aG_NEDSens) * dt; %Roll/Pitch Must be Zero

% Position propagation
rPred = rPrev + vPrev * dt;

% Quaternion propagation
qPred = (eye(4) + (dt/2) * (s2b(wBody)))*qPrev;
% State Vector propagation without process noise
xProp = [rPred; vPred; qPred];
end

