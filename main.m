function [] = main()
% Title:            Inertial Navigation
% Sub Title:        EKF Sensor Fusion (IMU/GPS)
% Script Writer:	Awais Arshad
% Association:      ASCL, KAIST
% Date:             June 29th, 2020
% Advisor:          Prof. Hyuchoong Bang
% email:            m.awais@kaist.ac.kr
% Code Ref:         https://openimu.readthedocs.io/en/latest/algorithms.html
% Useful Ref(s):    https://www.telesens.co/category/sensor-fusion/

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
global stdGpsPos stdGpsVel stdEuler
stdGyro = 0.00005; %rad/sec, Measurement Error
stdAcc = 0.00005; %m/sec^2, Measurement Error
stdDriftDotGyro = 0.001; % standard deviation of Gyro Drift
stdDriftDotAcc = 0.001; % Standard deviation of Acc Drift
% --> GyroBias = GyroOffSet + GyroDrift; GyroOffSet = zero
% --> AccBias = AccOffSet + AccDrift ; AccOffSet = zero

stdGpsPos = [0.1, 0.1, 0.3]; %Stand. Dev. of GPS Pos Measurement, meters 
stdGpsVel = [0.05, 0.05, 0.1]; %Std. dev. of GPS Vel Measurement, meters/sec
stdEuler = [0.01, 0.01, 0.01]; %Std. dev. of euler ang. meas. error, rad/sec


% --> Detailed specs are available in associated sensor measurement functions

% Simulation Parameters
N = 100;     % (dataSamples-1), 1st smaple used for initilization
M = 5;	% Number of Monte-Carlo runs
dt = 0.1;   % Sample Time

%% Extended Kalman Filter simulation
resXEst = zeros(16,N+1,M); % Monte-Carlo estimates
resXEstErr = zeros(16,N+1,M); % Monte-Carlo estimate errors
resPDiag = zeros(16,N+1); % diagonal term of estimation error covariance matrix
xTrue = zeros(19, N+1);

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
    xErrIni(7:10,1) = xTrueIni(7:10) + stdIni(3)*rand(4,1); %randQuat(); %
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
    for n = 1:1:N
        %%Prediction
        %Step 1a: State Prediction
        xPred = stateTransMdl(xPrev, dt, aMeasPrev, wMeasPrev);
        %Step 1b: Process Model
        [F, Q] = processMdl(xPrev, dt, aMeasPrev, wMeasPrev);

        %Step 2: Error Covariance Prediction (P -> Process estimation error covariance matrix)
        PPred = F*PPrev*F' + Q;
        
        %%Update
        %Step 3: Measurement Readout
        [zMeas, wMeas, aMeas, zTrue, R] = measSensorReading(n, gpsLLARef);
        H = observationMdl(xPred)
       
        %Step 4: Measurement Prediction
        h = measSensorPrediction(xPred);

        %Step 6: Kalman Gain
        S = H*PPred*H'+R;	%Innovation Matrix
        K = PPred*H'/S;     %Kalman Gain

        %Step 7: State Estimation/Update
        xEst(:,n) = xPred + K * (zMeas - h);
        if m==1
            xTrue(1:6,n) = zTrue(1:6,1);
            xTrue(7:10,n) = euler2quat(zTrue(7:9,1));
            xTrue(11:16,n) = zeros(6,1);
            xTrue(17:19,n) = zTrue(7:9,1);
        end

        %Step 8: Error Covariance Estimation/Update
        PEst = (eye(16) - K*H)*PPred;
        PDiag(:,n) = diag(PEst);
 
        xPrev = xEst(:,n);
        PPrev = PEst;
        wMeasPrev = wMeas;
        aMeasPrev = aMeas;
        if (m == 1)
            xTrue(11:16,n) = xEst(11:16,n); %Bug
        end
    end
    PDiagComp = [PDiagIni, PDiag];
    xEstComp = [xErrIni, xEst];
    if (m == 1)
        xTrueIni(17:19,1) = quat2euler(xTrueIni(7:10, 1)); %Data augmented for analysis only
        xTrue = [xTrueIni, xTrue(:,1:N)]; %For analysis only
    end
    resXEst(:,:,m) = xEstComp;
%    resXTrue(:,:,m) = xTrueComp;    %For analysis only
    resXEstErr(:,:,m) = xEstComp - xTrue(1:16,:);
    resPDiag(:,:,m) = PDiagComp;
end

%% get result statistics
xEstAVG = mean(resXEst,3); %Average of all monte carlo runs
xEstErrAVG = mean(resXEstErr,3); %Average of all monte carlo runs
x_RMSE = zeros(size(resXEstErr, 1),N+1); % root mean square error
PDiagAVG = mean(resPDiag,3);
% xTrue = mean(resXTrue,3); %Meaningless quantity. 4

for n = 1:1:N+1
    xEstAVG(17:19,n) = quat2euler(xEstAVG(7:10)); %Data augmented with euler angles for analysis purpose
    xEstErrAVG(17:19,n) = quat2euler(xEstErrAVG(7:10)); %Data augmented with euler angles for analysis purpose
    xTrue(17:19,n) = quat2euler(xTrue(7:10)); %Data augmented with euler angles for analysis purpose
    for m = 1:1:size(resXEstErr, 1)
        x_RMSE(m,n) = sqrt(mean(resXEstErr(m,n,:).^2,3));
    end
    x_RMSE(17:19,n) = quat2euler(x_RMSE(7:10));
end


%% plot results
time = (0:1:N)*dt;

% dispEstStates(NoOfSamples, sampleTime, stateNum)
NoOfSamples = N;
sampleTime = dt;

xPosIdx = 1:3;
xVelIdx = 4:6;
xQuatIdx = 7:10;
xWBIdx = 11:13;
xABIdx = 14:16;
xEulerIdx = 17:19;

dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xPosIdx); suptitle('True/Estimated Position, ');
dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xVelIdx); suptitle('True/Estimated Velocity, ');
dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xQuatIdx); suptitle('True/Estimated Quaternion, ');
dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xEulerIdx); suptitle('True/Estimated Euler, ');
dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xWBIdx); suptitle('True/Estimated Gyro Bias, ');
dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, xABIdx); suptitle('True/Estimated Accelerometer Bias, ');

dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xPosIdx);  suptitle('Standard Deviation of Position Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xVelIdx);  suptitle('Standard Deviation of Velocity Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xQuatIdx);  suptitle('Standard Deviation of Quaternion Errors, ');
% dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xEulerIdx);suptitle('Standard Deviation of Euler Errors, ');
% Covariance propagation from quaternion to euler is pending yet
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xWBIdx);  suptitle('Standard Deviation of Gyro Bias Errors, ');
dispEstErrors(x_RMSE, PDiagComp, NoOfSamples, sampleTime, xABIdx);  suptitle('Standard Deviation of Acclerometer Bias Errors, ');

dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xPosIdx); suptitle('Position Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xVelIdx); suptitle('Velocity Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xQuatIdx); suptitle('Quaternion Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xEulerIdx); suptitle('Euler Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
% Covariance propagation from quaternion to euler is pending yet
dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xWBIdx); suptitle('Acceleration Bias Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');
dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xABIdx); suptitle('Angular Velocity Error with 3-Sigma Bounds View and Absolute Accelerometer Bias');


dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xPosIdx); suptitle('Position Error with 3-Sigma Bounds View and Absolute Gyro Bias');
dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xVelIdx); suptitle('Velocity Error with 3-Sigma Bounds View and Absolute Gyro Bias');
dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xQuatIdx); suptitle('Quaternion Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xEulerIdx); suptitle('Euler Error with 3-Sigma Bounds View and Absolute Gyro Bias');
% Covariance propagation from quaternion to euler is pending yet
dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xWBIdx); suptitle('Acceleration Bias Error with 3-Sigma Bounds View and Absolute Gyro Bias');
dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErrAVG, PDiagComp, NoOfSamples, sampleTime, xABIdx); suptitle('Angular Velocity Error with 3-Sigma Bounds View and Absolute Gyro Bias');

dispEstStates3D(xTrue, xEstAVG, xPosIdx); suptitle('True/Estimated Position, ');
dispEstStates3D(xTrue, xEstAVG, xVelIdx); suptitle('True/Estimated Velocity, ');
% 


end


function dispEstStates(xTrue, xEstAVG, NoOfSamples, sampleTime, stateNum)
figure('Name', 'Time history of an estimation results');
time = (0:1:NoOfSamples)*sampleTime;
for n = 1:1:length(stateNum)
    subplot(length(stateNum),1,n); hold on;
    plot(time, xTrue(stateNum(n),:), 'linewidth', 2);
    plot(time, xEstAVG(stateNum(n),:), '--', 'linewidth', 2);
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



function dispEstStates3D(xTrue, xEstAVG, sVec3D)
% time = (0:1:NoOfSamples)*sampleTime;
figure('Name', 'estimation results comparison in 3D');
plot3(xTrue(sVec3D(1),:), xTrue(sVec3D(2),:), xTrue(sVec3D(3),:), '-o','MarkerSize',5); hold on;
plot3(xEstAVG(sVec3D(1),:), xEstAVG(sVec3D(2),:), xEstAVG(sVec3D(3),:), '-o','MarkerSize',5);
legend({'True', 'Estimated'}, 'fontsize', 12);
xlabel('X', 'fontsize', 12); ylabel('Y', 'fontsize', 12); zlabel('Z', 'fontsize', 12);
grid on;
end


function dispEstErr_3SigmaBound_AccBias(xEstAVG, xEstErr_AVG, PDiagComp, NoOfSamples, sampleTime, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Accelerometer Bias');
time = (0:1:NoOfSamples)*sampleTime;


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
    legend([e b1 biasAcc], {'Estimation Error', '3-Sigma Bounds', 'Accelerometer Bias'}, 'fontsize', 12);
    xlabel('Time (sec)', 'fontsize', 12);
    ylabel(strcat('X', num2str(stateNum(n), '%0.2d'), ' '), 'fontsize', 12); grid on;
end

end

function dispEstErr_3SigmaBound_GyroBias(xEstAVG, xEstErr_AVG, PDiagComp, NoOfSamples, sampleTime, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Gyro Bias');
time = (0:1:NoOfSamples)*sampleTime;


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

function dispEstErr_3SigmaBound_AccGyroBias(xEstAVG, xEstErr_AVG, PDiagComp, NoOfSamples, sampleTime, stateNum)
figure('Name','Estimate errors along with 3-sigma bounds and Accelerometer Bias');
time = (0:1:NoOfSamples)*sampleTime;


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

