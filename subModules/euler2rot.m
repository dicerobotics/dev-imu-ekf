function dcm = euler2rot( euler )
% Adapted from matlab function dcm2angle with rotation sequence ZYX
% Modified by:          Awais Arshad
% Modification Date:	July 9th, 2020

% Euler Format -> [roll, pitch, yaw];

angles = [euler(3) euler(2) euler(1)];

dcm = zeros(3,3);
cang = cos( angles );
sang = sin( angles );

dcm(1,1) = cang(2).*cang(:,1);
dcm(1,2) = cang(2).*sang(:,1);
dcm(1,3) = -sang(2);
dcm(2,1) = sang(3).*sang(2).*cang(1) - cang(3).*sang(1);
dcm(2,2) = sang(3).*sang(2).*sang(1) + cang(3).*cang(1);
dcm(2,3) = sang(3).*cang(2);
dcm(3,1) = cang(3).*sang(2).*cang(1) + sang(3).*sang(1);
dcm(3,2) = cang(3).*sang(2).*sang(1) - sang(3).*cang(1);
dcm(3,3) = cang(3).*cang(2);
end

