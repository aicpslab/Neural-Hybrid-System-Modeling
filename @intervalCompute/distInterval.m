function I = distInterval(I1,I2)
% compute distance between two intervals I1 and I2
%
% Syntax:
%    I = disInterval(I1,I2)
%
% Inputs:
%    I1 - interval I1
%    I2 - interval I2
%
% Outputs:
%    I - interval matrix, minimal and maximal values of each dimension

% Author:       Weiming Xiang
% Written:      03/16/2020
% Last update:  03/16/2020

%------------- BEGIN CODE --------------
[dim1,~] = size(I1);
[dim2,~] = size(I2);
if dim1 == dim2
    dim = dim1;
else
    error('Error. Two intervals should have same dimensions.')
end
for i = 1:1:dim
    d1 = abs(I2(i,2)-I1(i,1));
    d2 = abs(I1(i,2)-I2(i,1));
    if I1(i,2) < I2(i,1) || I2(i,2) < I1(i,1)
        I(i,1) = min(d1,d2);
        I(i,2) = max(d1,d2);
    else
        I(i,1) = 0;
        I(i,2) = max(d1,d2);
    end
end
        
    



%------------- END OF CODE --------------