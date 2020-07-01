function [] = main()
% Title:            Inertial Navigation
% Sub Title:        EKF Sensor Fusion (IMU/GPS)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html

% State Vector
% x = [r, v, q, wb, ab] State Vector, (1 x 16)
% r = [rx, ry, rz]      NED Position
% v = [vx, vy, vz]      NED Velocity
% q = [q0, q1, q2, q3]  Body Attitude, q0 is scalar
% wb = [wbx, wby, wbz]  Angular Rate Bias
% ab = [abx, aby, abz]  Accelerometer Bias

%% Initialization
clc; clear; close all;
% Include directories
addpath('./subModules');
addpath('./dataSet/oxts/data');

% NED frame origin definition
gpsLLARef = [49.011642946299 8.4160639632147 115.71427154541]; 

% Sensor specifications
global stdGyro stdAcc stdDriftDotGyro stdDriftDotAcc
stdGyro = 0.00005; %rad/sec, Measurement Error
stdAcc = 0.00005; %m/sec^2, Measurement Error
stdDriftDotGyro = 0.001; % standard deviation of Gyro Drift
stdDriftDotAcc = 0.001; % Standard deviation of Acc Drift
% --> GyroBias = GyroOffSet + GyroDrift; GyroOffSet = zero
% --> AccBias = AccOffSet + AccDrift ; AccOffSet = zero


% --> Detailed specs are available in associated sensor measurement functions

% Simulation Parameters
N = 100;     % (dataSamples-1), 1st smaple used for initilization
M = 10;	% Number of Monte-Carlo runs
dt = 0.1;   % Sample Time

%% Extended Kalman Filter simulation
resXEst = zeros(16,N+1,M); % Monte-Carlo estimates
resXEstErr = zeros(16,N+1,M); % Monte-Carlo estimate errors
resPDiag = zeros(16,N+1); % diagonal term of estimation error covariance matrix

% filtering
for m = 1:1:M
    % Filter Parameter Initialization
    % --> x = [r, v, q, wb, ab] State Vector, (16 x 1)
    [zMeas, wMeas, aMeas, zTrueIni, R] = measSensorReading(0, gpsLLARef);
    xTrueIni = [zTrueIni(1:6,1); euler2quat(zTrueIni(7:9)); zeros(6,1)];
    stdIni = [0.5, 0.25, 0, 1, 1]'; % std for Errpr in initial guess, 
%   stdIni = [Pos, Vel, quatScalar, QuatVector, biasW, biasA];
    xErrIni(1:3,1) = xTrueIni(1:3) + stdIni(1)*rand(3,1);
    xErrIni(4:6,1) = xTrueIni(4:6) + stdIni(2)*rand(3,1);
    xErrIni(7:10,1) = xTrueIni(7:10) + stdIni(3)*rand(4,1);
    xErrIni(11:13,1) = xTrueIni(11:13) + stdIni(4)*randn(3,1);
    xErrIni(14:16,1) = xTrueIni(14:16) + stdIni(5)*randn(3,1);
    xPrev = xErrIni;

    PPrev = [eye(3)*stdIni(1)^2, zeros(3,16-3);
              zeros(3,3), eye(3)*stdIni(2)^2, zeros(3, 16-6);
              zeros(4,6), eye(4)*stdIni(3)^2, zeros(4, 16-10);
              zeros(3,10), eye(3)*stdIni(4)^2, zeros(3,16-13);
              zeros(3, 13), eye(3)*stdIni(5)^2];

    PDiagIni = diag(PPrev);
    wMeasPrev = wMeas;
    aMeasPrev = aMeas;
    for k = 1:1:N
        %%Prediction
        %Step 1a: State Prediction
        xPred = stateTransMdl(xPrev, dt, aMeasPrev, wMeasPrev);
        %Step 1b: Process Model
        [F, Q] = processMdl(xPrev, dt, aMeasPrev, wMeasPrev);

        %Step 2: Error Covariance Prediction (P -> Process estimation error covariance matrix)
        PPred = F*PPrev*F' + Q;
        
        %%Update
        %Step 3: Measurement Readout
        [zMeas, wMeas, aMeas, zTrue, R] = measSensorReading(k, gpsLLARef);
         H = observationMdl(xPred);
       
        %Step 4: Measurement Prediction
        h = measSensorPrediction(xPred);

        %Step 6: Kalman Gain
        S = H*PPred*H'+R;	%Innovation Matrix
        K = PPred*H'/S;     %Kalman Gain

        %Step 7: State Estimation/Update
        xEst(:,k) = xPred + K * (zMeas - h);
        xTrue(:,k) = xPred + K * (zTrue - h); % Just for analysis

        %Step 8: Error Covariance Estimation/Update
        PEst = (eye(16) - K*H)*PPred;
        PDiag(:,k) = diag(PEst);
 
        xPrev = xEst(:,k);
        PPrev = PEst;
        wMeasPrev = wMeas;
        aMeasPrev = aMeas;
    end
    PDiagComp = [PDiagIni, PDiag];
    xEstComp = [xErrIni, xEst];
    xTrueComp = [xTrueIni, xTrue]; %For analysis only
    resXEst(:,:,m) = xEstComp;
    resXTrue(:,:,m) = xTrueComp;    %For analysis only
    resXEstErr(:,:,m) = xEstComp - xTrueComp;
    resPDiag(:,:,m) = PDiagComp;
end

%% get result statistics
xEstAVG = mean(resXEst,3); %Average of all monte carlo runs
xEstErr_AVG = mean(resXEstErr,3); %Average of all monte carlo runs
x_RMSE = zeros(size(resXEstErr, 1),N+1); % root mean square error
PDiagAVG = mean(resPDiag,3);
% xTrueAVG = mean(resXTrue,3); %Meaningless quantity. 4

for k = 1:1:N+1
    for m = 1:1:size(resXEstErr, 1)
        x_RMSE(m,k) = sqrt(mean(resXEstErr(m,k,:).^2,3));
    end
 end

%% plot results
time = (0:1:N)*dt;

% dispEstStates(NoOfSamples, sampleTime, stateNum)
NoOfSamples = N;
sampleTime = dt;
stateNum = [1:3];
xPosIdx = [1:3];
xVelIdx = [4:6];
xQuatIdx = [7:10];
xWBIdx = [11:13];
xABIdx = [14:16];

dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, xPosIdx); suptitle('True/Estimated Position, ');
dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, xVelIdx); suptitle('True/Estimated Velocity, ');
dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, xQuatIdx); suptitle('True/Estimated Quaternion, ');
dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, xWBIdx); suptitle('True/Estimated Gyro Bias, ');
dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, xABIdx); suptitle('True/Estimated Accelerometer Bias, ');

dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xPosIdx);  suptitle('Standard Deviation of Position Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xVelIdx);  suptitle('Standard Deviation of Velocity Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xQuatIdx);  suptitle('Standard Deviation of Quaternion Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xWBIdx);  suptitle('Standard Deviation of Gyro Bias Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xABIdx);  suptitle('Standard Deviation of Acclerometer Bias Errors, ');


dispEstStates3D(resXTrue, resXEst, xPosIdx); suptitle('True/Estimated Position, ');
dispEstStates3D(resXTrue, resXEst, xVelIdx); suptitle('True/Estimated Velocity, ');



end


function dispEstStates(resXTrue, resXEst, NoOfSamples, sampleTime, stateNum)
figure('Name', 'Time history of an estimation results');
time = (0:1:NoOfSamples)*sampleTime;
for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, resXTrue(stateNum(n),:, 1), 'linewidth', 2);
    plot(time, resXEst(stateNum(n),:,1), '--', 'linewidth', 2);
    legend({'True', 'Estimated'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d')), 'fontsize', 12); grid on;
end
end

function dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, stateNum)
figure('Name','Actual and estimated standard deviation for estimate errors');
time = (0:1:NoOfSamples)*sampleTime;


for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, x_RMSE(stateNum(n),:), 'linewidth', 2);
    plot(time, sqrt(PDiagComp(stateNum(n),:)), '--', 'linewidth', 2);
    legend({'RMSE Estimation', 'RMSE Process'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' Error Std'), 'fontsize', 12); grid on;
end

end



function dispEstStates3D(resXTrue, resXEst, sVec3D)
% time = (0:1:NoOfSamples)*sampleTime;
figure('Name', 'estimation results comparison in 3D');
plot3(resXTrue(sVec3D(1),:, 1), resXTrue(sVec3D(2),:, 1), resXTrue(sVec3D(3),:, 1), '-o','MarkerSize',5); hold on;
plot3(resXEst(sVec3D(1),:, 1), resXEst(sVec3D(2),:, 1), resXEst(sVec3D(3),:, 1), '-o','MarkerSize',5);
legend({'True', 'Estimated'}, 'fontsize', 12);
xlabel('X', 'fontsize', 12); ylabel('Y', 'fontsize', 12); zlabel('Z', 'fontsize', 12);
grid on;
end
