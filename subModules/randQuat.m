function out = randQuat()
% 
% while true
%     x = -1 + 2 * rand;
%     y = -1 + 2 * rand;
%     z = x*x + y*y;
%     if(z < 1) 
%         break;
%     end
% end
% 
% while true
%     u = -1 + 2 * rand;
%     v = -1 + 2 * rand;
%     w = u*u + v*v;
%     if(w < 1) 
%         break;
%     end
% end
% 
% s = sqrt((1-z) / w);
% out = [x, y, s*u, s*v]';

%%%%%%%
% u = -1 + 2 * rand;
% v = -1 + 2 * rand;
% w = -1 + 2 * rand;
% out = real([sqrt(1-u)*sin(2*pi*v), sqrt(1-u)*cos(2*pi*v), sqrt(u)*sin(2*pi*w), sqrt(u)*cos(2*pi*w)])';

% v = rand(3,1);
% v = v/norm(v);
% out = [rand; v];
theta = 2*pi*rand - pi;
v = rand(3,1);
vNorm = v/norm(v);

out = [cos(theta/2); vNorm*sin(theta/2)];
end
