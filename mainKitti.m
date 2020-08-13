clc; clear all; close all;% rng('Default');

% addpath('./dataSet/oxtsSync', './dataSet/oxtsSync/data');
% mainKittiOxts();
% rmpath('./dataSet/oxtsSync', './dataSet/oxtsSync/data');

addpath('./dataSet/oxtsUnSync', './dataSet/oxtsUnSync/data');
mainKittiOxts();
rmpath('./dataSet/oxtsUnSync', './dataSet/oxtsUnSync/data');

function [] = mainKittiOxts()
% Title:            Inertial Navigation
% Sub Title:        EKF Sensor Fusion (IMU/GPS)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             Aug 13th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Status:           Development in Progress
% Code Ref:         [1] {Kalman Filter} https://openimu.readthedocs.io/en/latest/algorithms.html
% Useful Ref(s):    [2] {Kalman Filter+Matlab Simulation} https://www.telesens.co/category/sensor-fusion/
%                   [3] {Inertial Nav + Noises} https://www.cl.cam.ac.uk/techreports/UCAM-CL-TR-696.pdf
%                   [4] {IMU Datasheet} https://www.oxts.com/app/uploads/2017/07/RT3000-brochure-170606.pdf
%                   [5] {Position in Ref-Frames} https://www.youtube.com/watch?v=GUvoVvXwoOQ&list=PLmiTuMrmGgEMmjb6IXGh8RMWaUV_cS-2v&index=1  
%                       (Time: 28:00 -> End)
%                   [6] {Velocity in Ref-Frames} https://www.youtube.com/watch?v=ZNVvYg1FOPk&list=PLmiTuMrmGgEMmjb6IXGh8RMWaUV_cS-2v&index=2 
%                       (Time: start -> 28:00)
%                   [7] {Quaternion Basics} http://www.tu-berlin.de/fileadmin/fg169/miscellaneous/Quaternions.pdf
%                   [8] {Kinematics of Moving Frames} https://ocw.mit.edu/courses/mechanical-engineering/2-017j-design-of-electromechanical-robotic-systems-fall-2009/course-text/MIT2_017JF09_ch09.pdf

%%
% State Vector
% x = [r, v, q, wb, ab] State Vector, (1 x 16)
% r = [rx, ry, rz]      NED Position
% v = [vx, vy, vz]      NED Velocity
% q = [q0, q1, q2, q3]  Body Attitude, q0 is scalar
% wb = [wbx, wb y, wbz]  Angular Rate Bias
% ab = [abx, aby, abz]  Accelerometer Bias

%% Initialization

% Include directories
addpath('./subModules');

%Read Data Time Stamps
fileName = 'timestamps.txt'; fileID = fopen(fileName, 'r');
formatSpec = '%4.0f-%2.0f-%2.0f %2.0f:%2.0f:%2.0f.%f';
timeStamps = cell2mat(textscan(fileID,formatSpec));
timeInSec = timeStamps(:, 4)*3600 + timeStamps(:, 5)*60 + timeStamps(:, 6) + timeStamps(:, 7) * 10^(-9); %Assumption: Data collected on the same date
timeSeries = timeInSec - ones(size(timeInSec))*timeInSec(1);
dtTrue = timeSeries(2:end) - timeSeries(1:end-1);

% NED frame origin definition
fileName = strcat(num2str(0, '%0.10d'), '.txt'); data = load(fileName);
gpsLLARef = [data(1), data(2), data(3)];
dtMean = mean(dtTrue);	%Sample Time

% Sensor specifications
global stdGyro stdAcc stdDriftDotGyro stdDriftDotAcc
global stdGpsPos stdGpsVel stdEuler
global gyroOffSet accOffSet

ARW_GyroDataSheet = 0.2;	%Angular Random Walk, Units: degree/rt-hr, OXTS Inertial+GNSS RT3000 v2 - (RT3003)
BS_GyroDataSheet = 2.0;     %Bias Stability, Units: degree/hr, OXTS Inertial+GNSS RT3000 v2 - (RT3003)
VRW_AccDataSheet = 0.005;	%Velocity Random Walk, Units: m/s/rt-hr, OXTS Inertial+GNSS RT3000 v2 - (RT3003)
BS_AccDataSheet = 2e-6;     %Accelerometer Bias Stability, Units: g,  OXTS Inertial+GNSS RT3000 v2 - (RT3003)

ARW_Gyro = (ARW_GyroDataSheet * (pi/180)) / (sqrt(3600));   % Units: rad/rt-sec
BS_Gyro = BS_GyroDataSheet * (pi/180) / (3600);             % Units: rad/sec
VRW_Acc = VRW_AccDataSheet * (1/sqrt(3600));                % Units: m/s/s
BS_Acc = BS_AccDataSheet * 9.80665;                         % Units: m/s

stdGyro = ARW_Gyro/sqrt(dtMean);    %Ref [3], Thermo-Mechanical White Noise -> Angle Random Walk
stdAcc = VRW_Acc/sqrt(dtMean);      %Ref [3], Thermo-Mechanical White Noise -> Angle Random Walk
stdDriftDotGyro = (2*pi/log(2))*(BS_Gyro^2/ARW_Gyro);	%Ref: https://openimu.readthedocs.io/en/latest/algorithms/STM_Bias.html
stdDriftDotAcc = (2*pi/log(2))* (BS_Acc^2/VRW_Acc); %Ref: https://openimu.readthedocs.io/en/latest/algorithms/STM_Bias.html

stdGyro = 1 * stdGyro;
stdDriftDotAcc = 1000 * stdDriftDotAcc; %Increase value see considerable effect in output graphs
stdDriftDotGyro = 10000 * stdDriftDotGyro; %Increase value see considerable effect in output graphs

% Note: Further literature survey suggested 

gyroOffSet = zeros(3,1);	% --> GyroBias = GyroOffSet + GyroDrift;
accOffSet = zeros(3,1);     % --> AccBias = AccOffSet + AccDrift;

stdGpsPos = 0.5 * ones(1,3);                           %CEP, Units: meters, Approximation: CEP=1-sigma,  OXTS Inertial+GNSS RT3000 v2 - (RT3003)
stdGpsVel = 0.05*(1000/(60*60))*ones(1,3);             %RMS, Units: m/s, Approximation: RMS=1-sigma, OXTS Inertial+GNSS RT3000 v2 - (RT3003)
stdEuler = ones(1,3) .* deg2rad([0.03, 0.03, 0.1]);    %1-sigma, Units: radian, OXTS Inertial+GNSS RT3000 v2 - (RT3003)

% Simulation Parameters
N = length(dtTrue); % (dataSamples-1), 1st smaple used for initilization
M = 1;              % Number of Monte-Carlo runs

%% Extended Kalman Filter simulation
resXEst = zeros(16,N+1,M);      % Monte-Carlo estimates
xMeas = zeros(19, N);
xEst = zeros(16, N);

% Filtering
for m = 1:1:M
    % Filter Parameter Initialization
    % --> x = [r, v, q, wb, ab] State Vector, (16 x 1)
    [zMeas, wMeas, aMeas, ~] = measSensorReading(0, gpsLLARef, dtMean);

    xMeasInit = [zMeas(1:6,1); angle2quat(zMeas(9), zMeas(8), zMeas(7), 'ZYX')'; gyroOffSet; accOffSet];
    xTrueInit = [zMeas(1:6,1); angle2quat(zMeas(9), zMeas(8), zMeas(7), 'ZYX')';  gyroOffSet; accOffSet];

    stdIni = [mean(stdGpsPos), mean(stdGpsVel), 0, 0, 0]'; % std for Error in initial guess, %stdIni -> [Pos, Vel, euler, biasW, biasA];
    xInit(1:3,1) = xTrueInit(1:3) + stdIni(1)*randn(3,1);
    xInit(4:6,1) = xTrueInit(4:6) + stdIni(2)*randn(3,1);
    xInit(7:10,1) = angle2quat(zMeas(9), zMeas(8), zMeas(7), 'ZYX')';
    xInit(11:13,1) = xTrueInit(11:13);% + stdIni(4)*randn(3,1);%zeros(3,1);%
    xInit(14:16,1) = xTrueInit(14:16);% + stdIni(5)*randn(3,1);%zeros(3,1);%

    PInit = [eye(3)*stdIni(1)^2, zeros(3,16-3);
              zeros(3,3), eye(3)*stdIni(2)^2, zeros(3, 16-6);
              zeros(4,6), eye(4)*stdIni(3)^2, zeros(4, 16-10);
              zeros(3,10), eye(3)*stdIni(4)^2, zeros(3,16-13);
              zeros(3, 13), eye(3)*stdIni(5)^2];    
    
    xPrev = xInit;
    PPrev = PInit;
    wMeasPrev = wMeas;
    aMeasPrev = aMeas;
   
for n = 1:1:N
        %%Prediction 
        %Step 1a: State Prediction
        xPred = stateTransMdl(xPrev, dtTrue(n), aMeasPrev, wMeasPrev);
        xPred = repairQuaternion(xPred);
        %Step 1b: Process Model
        [F, Q] = processMdl(xPrev, dtTrue(n), aMeasPrev, wMeasPrev);

        %Step 2: Error Covariance Prediction (P -> Process estimation error covariance matrix)
        PPred = F*PPrev*F' + Q; %Q matrix account for bias estimation
        
        %%Update
        %Step 3: Measurement Readout
        [zMeas, wMeas, aMeas, R] = measSensorReading(n, gpsLLARef, dtMean);
        H = obsMdl(xPred);
        %Step 4: Measurement Prediction
        h = measSensorPrediction(xPred);

        %Step 6: Kalman Gain
        S = H*PPred*H'+R;	%Innovation Matrix
        K = PPred*H'/S;     %Kalman Gain

        %Step 7: State Estimation/Update
        inovError = (zMeas - h); %Innovation Error
        xEst(:,n) = xPred + K * inovError; %zMeas Contains euler angles from NED to Body frame
        xEst(:,n) = repairQuaternion(xEst(:,n));% disp(xEst(11:16,n));
        if (m==1)
            xMeas(1:6,n) = zMeas(1:6,1); 
            xMeas(7:10,n) = angle2quat(zMeas(9,1), zMeas(8,1), zMeas(7,1), 'ZYX')'; 
            xMeas(11:16,n) = xMeasInit(11:16,1);%[gyroOffSet; accOffSet];%Bug: We don't have these values in dataset
            xMeas(17:19,n) = zMeas(7:9,1);
        end

        %Step 8: Error Covariance Estimation/Update
        PEst = (eye(16) - K*H)*PPred;
 
        xPrev = xEst(:,n);
        PPrev = PEst;
        wMeasPrev = wMeas;
        aMeasPrev = aMeas;
end
    xEstComp = [xInit, xEst];
    if (m == 1)
        [yaw, pitch, roll] = quat2angle(xMeasInit(7:10, 1)', 'ZYX');
        xMeasInit(17:19,1) = [roll, pitch, yaw]';
        xMeasComp = [xMeasInit, xMeas]; %For analysis only
        
    end
    resXEst(:,:,m) = xEstComp;
end

%%
xEstAVG = mean(resXEst,3); %Average of all monte carlo runs

for n = 1:1:N+1
    [yaw, pitch, roll] = quat2angle(xEstAVG(7:10,n)', 'ZYX'); xEstAVG(17:19,n) = [roll, pitch, yaw]';
end


%% plot results
time = timeSeries;
xPosIdx = 1:3; xVelIdx = 4:6; xQuatIdx = 7:10; xWBIdx = 11:13; xABIdx = 14:16; xEulerIdx = 17:19;

dispEstStates(xMeasComp, xEstAVG, xPosIdx, time); suptitle('Measured/Estimated Position');
dispEstStates(xMeasComp, xEstAVG, xVelIdx, time); suptitle('Measured/Estimated Velocity');
dispEstStates(xMeasComp, xEstAVG, xABIdx, time); suptitle('Estimated aBias, Ignore Measured Values');
dispEstStates(xMeasComp, xEstAVG, xWBIdx, time); suptitle('Estimated wBias, Ignore Measured Values');
dispEstStates(xMeasComp, xEstAVG, xQuatIdx, time); suptitle('Measured/Estimated Quaternion');
dispEstStates(xMeasComp, xEstAVG, xEulerIdx, time); suptitle('Measured/Estimated Euler');

dispEstStates3D(xMeasComp, xEstAVG, xPosIdx); suptitle('Measured/Estimated Position, ');
dispEstStates3D(xMeasComp, xEstAVG, xVelIdx); suptitle('Measured/Estimated Velocity, ');

rmpath('./subModules');
end


function dispEstStates(xMeas, xEstAVG, stateNum, time)
figure('Name', 'Time history of an estimation results');
for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, xMeas(stateNum(n),:), 'r-.', 'linewidth', 2);
    plot(time, xEstAVG(stateNum(n),:), 'b--', 'linewidth', 2);
    legend({'Measured', 'Estimated'}, 'fontsize', 12);
    
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d')), 'fontsize', 12); grid on;
end
end


function dispEstStates3D(xMeasComp, xEstAVG, sVec3D)
figure('Name', 'estimation results comparison in 3D');
plot3(xMeasComp(sVec3D(1),:), xMeasComp(sVec3D(2),:), xMeasComp(sVec3D(3),:), 'r-.', 'linewidth', 2); hold on;
plot3(xEstAVG(sVec3D(1),:), xEstAVG(sVec3D(2),:), xEstAVG(sVec3D(3),:), 'b--', 'linewidth', 2);

legend({'Measured', 'Estimated'}, 'fontsize', 12);
xlabel('X', 'fontsize', 12); ylabel('Y', 'fontsize', 12); zlabel('Z', 'fontsize', 12);
grid on;
end
