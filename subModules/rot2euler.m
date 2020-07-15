function euler = rot2euler( dcm)
% Adapted from matlab function dcm2angle with rotation sequence ZYX
% Modified by:    	Awais Arshad
% Modification Date:  July 9th, 2020
dcm = dcm';
[yaw,pitch,roll] = threeaxisrot( dcm(1,2), dcm(1,1), -dcm(1,3), ...
                                   dcm(2,3), dcm(3,3));
                              
euler = [roll, pitch, yaw]';

function [r1,r2,r3] = threeaxisrot(r11, r12, r21, r31, r32)
    % find angles for rotations about X, Y, and Z axes
    r1 = atan2( r11, r12 );
    r21(r21 < -1) = -1;
    r21(r21 > 1) = 1;
    r2 = asin( r21 );
    r3 = atan2( r31, r32 );
end

end
