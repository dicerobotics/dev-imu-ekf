function [] = main()
% Title:            Inertial Navigation
% Sub Title:        EKF Sensor Fusion (IMU/GPS)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
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

%%
% State Vector
% x = [r, v, q, wb, ab] State Vector, (1 x 16)
% r = [rx, ry, rz]      NED Position
% v = [vx, vy, vz]      NED Velocity
% q = [q0, q1, q2, q3]  Body Attitude, q0 is scalar
% wb = [wbx, wb y, wbz]  Angular Rate Bias
% ab = [abx, aby, abz]  Accelerometer Bias
%% Initialization
clc; clear; close all;% rng('Default');
% Include directories
addpath('./subModules');
addpath('./dataSet/oxts');
addpath('./dataSet/oxts/data');

%Read Data Time Stamps
fileName = 'timestamps.txt'; fileID = fopen(fileName, 'r');
formatSpec = '%4.0f-%2.0f-%2.0f %2.0f:%2.0f:%2.0f.%f';
timeStamps = cell2mat(textscan(fileID,formatSpec));
timeInSec = timeStamps(:, 4)*3600 + timeStamps(:, 5)*60 + timeStamps(:, 6) + timeStamps(:, 7) * 10^(-9); %Assumption: Data collected on the same date
timeSeries = timeInSec - ones(size(timeInSec))*timeInSec(1);
dtTrue = timeSeries(2:end) - timeSeries(1:end-1);

% NED frame origin definition
gpsLLARef = [49.008644826538 8.3981039999565 112.99059295654]; %lat,long,alt in file 0000000000.txt
dtMean = mean(dtTrue);   % Sample Time

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
BS_Acc = BS_AccDataSheet * 9.8;                             % Units: m/s

stdGyro = ARW_Gyro/sqrt(dtMean);    %Ref [3], Thermo-Mechanical White Noise -> Angle Random Walk
stdAcc = VRW_Acc/sqrt(dtMean);      %Ref [3], Thermo-Mechanical White Noise -> Angle Random Walk
stdDriftDotGyro = (2*pi/log(2))*(BS_Gyro^2/ARW_Gyro);	%Ref: https://openimu.readthedocs.io/en/latest/algorithms/STM_Bias.html
stdDriftDotAcc = (2*pi/log(2))* (BS_Acc^2/VRW_Acc); %Bug Suspected   %Ref: https://openimu.readthedocs.io/en/latest/algorithms/STM_Bias.html
% Note: Further literature survey suggested 

gyroOffSet = zeros(3,1);	% --> GyroBias = GyroOffSet + GyroDrift;
accOffSet = zeros(3,1);     % --> AccBias = AccOffSet + AccDrift ;

stdGpsPos = 0.5 * ones(1,3);                           %CEP, Units: meters, Approximation: CEP=1-sigma,  OXTS Inertial+GNSS RT3000 v2 - (RT3003)
stdGpsVel = 0.05*(1000/(60*60))*ones(1,3);             %RMS, Units: m/s, Approximation: RMS=1-sigma, OXTS Inertial+GNSS RT3000 v2 - (RT3003)
stdEuler = ones(1,3) .* deg2rad([0.03, 0.03, 0.1]);    %1-sigma, Units: radian, OXTS Inertial+GNSS RT3000 v2 - (RT3003)

% Simulation Parameters
N = length(dtTrue); % (dataSamples-1), 1st smaple used for initilization
M = 1;              % Number of Monte-Carlo runs

%% Extended Kalman Filter simulation
resXEst = zeros(16,N+1,M);      % Monte-Carlo estimates
resXEstErr = zeros(16,N+1,M);   % Monte-Carlo estimate errors
resPDiag = zeros(16,N+1);       % Diagonal term of estimation error covariance matrix
xMeas = zeros(19, N);

% Filtering
for m = 1:1:M
    % Filter Parameter Initialization
    % --> x = [r, v, q, wb, ab] State Vector, (16 x 1)
    [zMeas, wMeas, aMeas, ~] = measSensorReading(0, gpsLLARef);
    xMeasInit = [zMeas(1:6,1); euler2quat(zMeas(7:9)); gyroOffSet; accOffSet];
    stdIni = 1*[0.2, 0.01, 0.1, 0.001, 0.001]'; % std for Error in initial guess, %stdIni -> [Pos, Vel, euler, biasW, biasA];
    xInit(1:3,1) = xMeasInit(1:3) + stdIni(1)*randn(3,1);
    xInit(4:6,1) = xMeasInit(4:6) + stdIni(2)*randn(3,1);
    xInit(7:10,1) = euler2quat(quat2euler(xMeasInit(7:10))+ stdIni(3)* randn(3,1));%randQuat(); % xTrueInit(7:10) + euler2quat(stdIni(3)* randn(3,1)); %
    xInit(7:10,1) = repairQuaternion(xInit(7:10,1));
    xInit(11:13,1) = xMeasInit(11:13) + stdIni(4)*randn(3,1);
    xInit(14:16,1) = xMeasInit(14:16) + stdIni(5)*randn(3,1);

    PInit = [eye(3)*stdIni(1)^2, zeros(3,16-3);
              zeros(3,3), eye(3)*stdIni(2)^2, zeros(3, 16-6);
              zeros(4,6), eye(4)*stdIni(3)^2, zeros(4, 16-10);
              zeros(3,10), eye(3)*stdIni(4)^2, zeros(3,16-13);
              zeros(3, 13), eye(3)*stdIni(5)^2];    
    
    PDiagInit = diag(PInit);
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
        PPred = F*PPrev*F' + Q;
        
        %%Update
        %Step 3: Measurement Readout
        [zMeas, wMeas, aMeas, R] = measSensorReading(n, gpsLLARef);
        H = symObservationMdl(xPred);
        %Step 4: Measurement Prediction
        h = measSensorPrediction(xPred);

        %Step 6: Kalman Gain
        S = H*PPred*H'+R;	%Innovation Matrix
        K = PPred*H'/S;     %Kalman Gain

        %Step 7: State Estimation/Update
        xEst(:,n) = xPred + K * (zMeas - h);
        xEst(:,n) = repairQuaternion(xEst(:,n));% disp(xEst(11:16,n));
        if (m==1)
            xMeas(1:6,n) = zMeas(1:6,1);
            xMeas(7:10,n) = euler2quat(zMeas(7:9,1));
            xMeas(11:16,n) = zeros(6,1); %Bug
            xMeas(17:19,n) = zMeas(7:9,1);
        end

        %Step 8: Error Covariance Estimation/Update
        PEst = (eye(16) - K*H)*PPred;
        PDiag(:,n) = diag(PEst);
 
        xPrev = xEst(:,n);
        PPrev = PEst;
        wMeasPrev = wMeas;
        aMeasPrev = aMeas;
        if (m == 1)
            xMeas(11:16,n) = xEst(11:16,n); %Bug, We don't have these values in dataset
        end
    end
    PDiagComp = [PDiagInit, PDiag];
    xEstComp = [xInit, xEst];
    if (m == 1)
        xMeasInit(17:19,1) = quat2euler(xMeasInit(7:10, 1)); %Data augmented for analysis only
        xMeasComp = [xMeasInit, xMeas]; %For analysis only
    end
    resXEst(:,:,m) = xEstComp;
    resXEstErr(:,:,m) = xEstComp - xMeasComp(1:16,:);
    resPDiag(:,:,m) = PDiagComp;
end

%% get result statistics
xEstAVG = mean(resXEst,3); %Average of all monte carlo runs
xEstErrAVG = mean(resXEstErr,3); %Average of all monte carlo runs
x_RMSE = zeros(size(resXEstErr, 1),N+1); % root mean square error
PDiagAVG = mean(resPDiag,3);

for n = 1:1:N+1
    xEstAVG(17:19,n) = quat2euler(xEstAVG(7:10,n)); %Data augmented with euler angles for analysis purpose
    xEstErrAVG(17:19,n) = quat2euler(xEstErrAVG(7:10,n)); %Data augmented with euler angles for analysis purpose
    for m = 1:1:size(resXEstErr, 1)
        x_RMSE(m,n) = sqrt(mean(resXEstErr(m,n,:).^2,3));
    end
    x_RMSE(17:19,n) = quat2euler(x_RMSE(7:10,n));
end


%% plot results
time = timeSeries;%(0:1:N)*dtMean;
NoOfSamples = N; %Excluding initial sample at t0
sampleTime = dtMean;

xPosIdx = 1:3;
xVelIdx = 4:6;
xQuatIdx = 7:10;
xWBIdx = 11:13;
xABIdx = 14:16;
xEulerIdx = 17:19;

dispEstStates(xMeasComp, xEstAVG, xPosIdx, time); suptitle('Measured/Estimated Position, ');
dispEstStates(xMeasComp, xEstAVG, xVelIdx, time); suptitle('Measured/Estimated Velocity, ');
dispEstStates(xMeasComp, xEstAVG, xQuatIdx, time); suptitle('Measured/Estimated Quaternion, ');
dispEstStates(xMeasComp, xEstAVG, xEulerIdx, time); suptitle('Measured/Estimated Euler, ');


% dispEstErrors(x_RMSE, PDiagComp, xPosIdx, time);  suptitle('Standard Deviation of Position Errors, ');
% dispEstErrors(x_RMSE, PDiagComp, xVelIdx, time);  suptitle('Standard Deviation of Velocity Errors, ');
% dispEstErrors(x_RMSE, PDiagComp, xQuatIdx, time);  suptitle('Standard Deviation of Quaternion Errors, ');
% % dispEstErrors(x_RMSE, PDiagComp, xEulerIdx, time);suptitle('Standard Deviation of Euler Errors, ');
% % Covariance propagation from quaternion to euler is pending yet
% dispEstErrors(x_RMSE, PDiagComp, xWBIdx, time);  suptitle('Standard Deviation of Gyro Bias Errors, ');
% dispEstErrors(x_RMSE, PDiagComp, xABIdx, time);  suptitle('Standard Deviation of Acclerometer Bias Errors, ');

dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xPosIdx); suptitle('Position Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xVelIdx); suptitle('Velocity Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xQuatIdx); suptitle('Quaternion Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% % dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xEulerIdx); suptitle('Euler Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% % Covariance propagation from quaternion to euler is pending yet
% dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xWBIdx); suptitle('Acceleration Bias Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, time, xABIdx); suptitle('Angular Velocity Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xPosIdx); suptitle('Position Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xVelIdx); suptitle('Velocity Error with 3-Sigma Bounds View and Absolute Gyro Bias');
dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xQuatIdx); suptitle('Quaternion Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% % dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xEulerIdx); suptitle('Euler Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% % Covariance propagation from quaternion to euler is pending yet
% dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xWBIdx); suptitle('Acceleration Bias Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xABIdx); suptitle('Angular Velocity Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% 
% dispEstStates3D(xMeasComp, xEstAVG, xPosIdx); suptitle('True/Estimated Position, ');
% dispEstStates3D(xMeasComp, xEstAVG, xVelIdx); suptitle('True/Estimated Velocity, ');
% 
% dispEstErr_3SigmaBound_AccGyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xPosIdx); suptitle('Position Error with 3-Sigma Bounds View and Biases');
% dispEstErr_3SigmaBound_AccGyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xVelIdx); suptitle('Velocity Error with 3-Sigma Bounds View and Biases');
% dispEstErr_3SigmaBound_AccGyroBias(xEstAVG, xEstErrAVG, PDiagComp, time, xQuatIdx); suptitle('Quaternion Error with 3-Sigma Bounds View and Biases');
% 
% misc(pTrueNED, time, 1:3); suptitle('Linear Position by dead reckoning');

end

function misc(xTrueComp, time, stateNum)
figure('Name', 'Time history of Acc/Gyro');
% time = (1:1:NoOfSamples)*sampleTime;
for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time(2:end), xTrueComp(stateNum(n),:), 'linewidth', 2);
    legend({'True'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel('Misc.', 'fontsize', 12); grid on;
end
end


function dispEstStates(xTrue, xEstAVG, stateNum, time)
figure('Name', 'Time history of an estimation results');
for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, xTrue(stateNum(n),:), 'linewidth', 2);
    plot(time, xEstAVG(stateNum(n),:), '--', 'linewidth', 2);
    legend({'Measured', 'Estimated'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d')), 'fontsize', 12); grid on;
end
end

function dispEstErrors(x_RMSE, PDiagComp, stateNum, time)
figure('Name','Actual and estimated standard deviation for estimate errors');
% time = (0:1:NoOfSamples)*sampleTime;


for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, x_RMSE(stateNum(n),:), 'linewidth', 2);
    plot(time, sqrt(PDiagComp(stateNum(n),:)), '--', 'linewidth', 2);
    
    legend({'RMSE Estimation', 'RMSE Process'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' Error Std'), 'fontsize', 12); grid on;
end

end



function dispEstStates3D(xTrue, xEstAVG, sVec3D)
% time = (0:1:NoOfSamples)*sampleTime;
figure('Name', 'estimation results comparison in 3D');
plot3(xTrue(sVec3D(1),:), xTrue(sVec3D(2),:), xTrue(sVec3D(3),:), '-o','MarkerSize',5); hold on;
plot3(xEstAVG(sVec3D(1),:), xEstAVG(sVec3D(2),:), xEstAVG(sVec3D(3),:), '-o','MarkerSize',5);
legend({'True', 'Estimated'}, 'fontsize', 12);
xlabel('X', 'fontsize', 12); ylabel('Y', 'fontsize', 12); zlabel('Z', 'fontsize', 12);
grid on;
end


function dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErr_AVG, PDiagComp, time, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Accelerometer Bias');
% time = (0:1:NoOfSamples)*sampleTime;


for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    e = plot(time, xEstErr_AVG(stateNum(n),:), 'k-',  'linewidth', 2);
    b1 = plot(time, 3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);
    plot(time, -3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);

    if (n >= 1) && (n <= 3)         %Position
        biasAcc = plot(time, xEstAVG(13+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 4) && (n <= 6)     %Velocity
        biasAcc = plot(time, xEstAVG(10+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 7) && (n <= 10)    %Quaternion
        biasAcc = plot(time, xEstAVG(7+n,:), 'b--', 'linewidth', 2); %Bug. No need to plot Acc Bias here
    elseif (n >= 11) && (n <= 13)   %Accelerometer Bias
        biasAcc = plot(time, xEstAVG(3+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 14) && (n <= 16)   %Gyro Bias
        biasAcc = plot(time, xEstAVG(0+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 17) && (n <= 19)   %Euler Angles
        biasAcc = plot(time, xEstAVG(-3+n,:), 'b--', 'linewidth', 2);
    end
    legend([e b1 biasAcc], {'Estimation Error', '3-Sigma Bounds', 'Accelerometer Bias'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' '), 'fontsize', 12); grid on;
end

end

function dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErr_AVG, PDiagComp, time, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Gyro Bias');
% time = (0:1:NoOfSamples)*sampleTime;


for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    e = plot(time, xEstErr_AVG(stateNum(n),:), 'k-',  'linewidth', 2);
    b1 = plot(time, 3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);
    plot(time, -3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);

    if (n >= 1) && (n <= 3)         %Position
        biasGyro = plot(time, xEstAVG(10+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 4) && (n <= 6)     %Velocity
        biasGyro = plot(time, xEstAVG(7+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 7) && (n <= 10)    %Quaternion
        biasGyro = plot(time, xEstAVG(4+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 11) && (n <= 13)   %Accelerometer Bias
        biasGyro = plot(time, xEstAVG(0+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 14) && (n <= 16)   %Gyro Bias
        biasGyro = plot(time, xEstAVG(-3+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 17) && (n <= 19)   %Euler Angles
        biasGyro = plot(time, xEstAVG(-6+n,:), 'b--', 'linewidth', 2);
    end
    legend([e b1 biasGyro], {'Estimation Error', '3-Sigma Bounds', 'Gyro Bias'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' '), 'fontsize', 12); grid on;
end

end

function dispEstErr_3SigmaBound_AccGyroBias(xEstAVG, xEstErr_AVG, PDiagComp, time, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Accelerometer Bias');
% time = (0:1:NoOfSamples)*sampleTime;


for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    e = plot(time, xEstErr_AVG(stateNum(n),:), 'k-',  'linewidth', 2);
    b1 = plot(time, 3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);
    plot(time, -3*sqrt(PDiagComp(stateNum(n),:)), 'r:', 'linewidth', 2);

    if (n >= 1) && (n <= 3)         %Position
        biasAcc = plot(time, xEstAVG(13+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 4) && (n <= 6)     %Velocity
        biasAcc = plot(time, xEstAVG(10+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 7) && (n <= 10)    %Quaternion
        biasAcc = plot(time, xEstAVG(7+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 11) && (n <= 13)   %Accelerometer Bias
        biasAcc = plot(time, xEstAVG(3+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 14) && (n <= 16)   %Gyro Bias
        biasAcc = plot(time, xEstAVG(0+n,:), 'b--', 'linewidth', 2);
    elseif (n >= 17) && (n <= 19)   %Euler Angles
        biasAcc = plot(time, xEstAVG(-3+n,:), 'b--', 'linewidth', 2);
    end
    
        if (n >= 1) && (n <= 3)         %Position
        biasGyro = plot(time, xEstAVG(10+n,:), 'm--', 'linewidth', 2);
    elseif (n >= 4) && (n <= 6)     %Velocity
        biasGyro = plot(time, xEstAVG(7+n,:), 'm--', 'linewidth', 2);
    elseif (n >= 7) && (n <= 10)    %Quaternion
        biasGyro = plot(time, xEstAVG(4+n,:), 'm--', 'linewidth', 2);
    elseif (n >= 11) && (n <= 13)   %Accelerometer Bias
        biasGyro = plot(time, xEstAVG(0+n,:), 'm--', 'linewidth', 2);
    elseif (n >= 14) && (n <= 16)   %Gyro Bias
        biasGyro = plot(time, xEstAVG(-3+n,:), 'm--', 'linewidth', 2);
    elseif (n >= 17) && (n <= 19)   %Euler Angles
        biasGyro = plot(time, xEstAVG(-6+n,:), 'm--', 'linewidth', 2);
        end
    
        
    legend([e b1 biasAcc biasGyro], {'Estimation Error', '3-Sigma Bounds', 'Accelerometer Bias', 'Gyro Bias'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' '), 'fontsize', 12); grid on;
end
end

