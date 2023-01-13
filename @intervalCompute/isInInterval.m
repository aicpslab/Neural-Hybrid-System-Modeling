function res = isInInterval(x,I)
% check  a point if inside an interval
%
% Syntax:
%    res = isInInterval(x,I)
%
% Inputs:
%    x - point
%    I - interval matrix
%
% Outputs:
%    res - indicate the interval is inside/outside/intersect
%    res = 0: I outside I_target
%    res = 1: I inside I_target

% Author:       Weiming Xiang
% Written:      02/25/2019
% Last update:  02/25/2019

%------------- BEGIN CODE --------------

[dim,~] = size(I);
res = 1;
for i = 1:1:dim
    if x(i) < I(i,1) || x(i) > I(i,2)
        res = 0;
        break
    end
end


    

%------------- END OF CODE --------------