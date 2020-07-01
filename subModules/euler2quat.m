function [quat] = euler2quat(euler)
%Ref: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
roll = euler(1);
pitch = euler(2);
yaw = euler(3);

cy = cos(yaw * 0.5);
sy = sin(yaw * 0.5);
cp = cos(pitch * 0.5);
sp = sin(pitch * 0.5);
cr = cos(roll * 0.5);
sr = sin(roll * 0.5);

q = zeros(4,1);
q(1) = cy * cp * cr + sy * sp * sr;
q(2) = cy * cp * sr - sy * sp * cr;
q(3) = sy * cp * sr + cy * sp * cr;
q(4) = sy * cp * cr - cy * sp * sr;
quat = q;

end