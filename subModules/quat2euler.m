function euler = quat2euler(quat)
q0 = quat(1); q1 = quat(2); q2 = quat(3); q3 = quat(4);
roll = atan2(2*(q2*q3+q0*q1), q0^2-q1^2-q2^2+q3^2);
pitch = -asin(2*(q1*q3-q0*q2));
yaw = atan2(2*(q1*q2+q0*q3), (q0^2+q1^2-q2^2-q3^2));
euler = [roll, pitch, yaw]';
end