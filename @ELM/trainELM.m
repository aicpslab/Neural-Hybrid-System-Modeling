function elm = trainELM(ELM1,Input,Output)
% ELMTRAIN Create and Train a Extreme Learning Machine
% Syntax
% [IW,B,LW,TF,TYPE] = elmtrain(ELM,Input,Output)
% Description
% Input  - Input Matrix of Training Set  (R*Q)
% Output  - Output Matrix of Training Set (S*Q)
% N   - Number of Hidden Neurons (default = Q)
% TF  - Transfer Function:
%       'sig' for Sigmoidal function (default)
%       'sin' for Sine function
%       'hardlim' for Hardlim function
% TYPE - Regression (0,default) or Classification (1)
% Output
% IW  - Input Weight Matrix (N*R)
% B   - Bias Matrix  (N*1)
% LW  - Layer Weight Matrix (N*S)
% Example
% Regression:
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',0)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)
% Classification
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',1)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)
% See also ELMPREDI
% Yu Lei,11-7-2010
% Copyright www.matlabsky.com
% $Revision:1.0 $
% if nargin < 5
%      TYPE = 0;
% if nargin < 4
%     TF = 'sig';
%      if nargin < 3
%         N = size(P,2);
%          if nargin < 2
%             error('ELM:Arguments','Not enough input arguments.');
%          end
%      end
% end   
% end

P=Input;
T=Output;
if size(P,2) ~= size(T,2)
    error('ELM:Arguments','The columns of P and T must be same.');
end
[R,~] = size(P);
% if TYPE  == 1
%     T  = ind2vec(T);
% end
[S,Q] = size(T);
if nargin < 2
     error('ELM:Arguments','Not enough input arguments.');
end
if size(ELM1.weight{1},2)~=R
    error('ELM:input dimension are not equal with ELM .');
end
if size(ELM1.weight{2},1)~=S
    error('ELM:output dimension are not equal with ELM .');
end
tic
BiasMatrix = repmat(ELM1.bias{1},1,Q);
% Calculate the Layer Output Matrix H
tempH = ELM1.weight{1} * P + BiasMatrix;
switch ELM1.activeFcn{1}
    case 'sig'
        H = 1 ./ (1 + exp(-tempH));
    case 'sin'
        H = sin(tempH);
    case 'hardlim'
        H = hardlim(tempH);
    case 'tansig'
        H = tansig(tempH);
    case 'ReLu'
        H = max(tempH,0);
end
% Calculate the Output Weight Matrix
LW = pinv(H') * T';
fprintf('time of training ELM is %d', toc)
fprintf(' seconds \n')
ELM1.weight{2}=LW';
X=LW'*H-T;
elm=ELM(ELM1.weight,ELM1.bias,ELM1.activeFcn);
elm.trainingError = partitions.MeanSquare(X);
end