function ELM1 = trainELMLipridge(ELM1,Input,Output)
%(P,T,N,TF,TYPE,d,delta)
% ELMTRAIN Create and Train a Extreme Learning Machine
% Syntax
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,N,TF,TYPE£¬tol)
% Description
% Input
% P   - Input Matrix of Training Set  (R*Q)
% T   - Output Matrix of Training Set (S*Q)
% N   - Number of Hidden Neurons (default = Q)
% TF  - Transfer Function:
%       'sig' for Sigmoidal function (default)
%       'sin' for Sine function
%       'hardlim' for Hardlim function
% TYPE - Regression (0,default) or Classification (1)
% Output
% IW  - Input Weight Matrix (N*R)
% B   - Bias Matrix  (N*1)
% LW  - Layer Weight Matrix (S*N)
% Example
% Regression:
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',0,tol)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)
% Classification
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',1)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)

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
IW=ELM1.weight{1};
%% Initial value of LW
%LW = pinv(H') * T';


delta=0;
for i = 1:size(Input,2)-1
    temp_delta=norm(Input(:,i)-Input(:,i+1));
    if (delta<temp_delta)
        delta=temp_delta;
    end 
end
% Y = elmpredict(P,ELMNetwork);
% figure('NumberTitle', 'off', 'Name', 'Initial Trained Erros')
% plot(P,Y-T);
% xlabel('P');
% ylabel('Initial Trained Erros');
%% Initial Network
% ELMNetwork1=ELMNetwork;
%% SPSA Algorithm to evalue the LW s.t min dist.||nn-ELMNetwork||
% LW=zeros(N,S);
% for n = 1:size(ELMNetwork.weight{2},1)
% %tic;
% %Lipshcitz & Ridge Regression 
%     %LW=SPSA(P,NN_ouputset,ELMNetwork,val,n,d);
     lambda = delta*norm(IW)*1*sqrt(2)*0.25;

%     LW(:,n) = pinv(H*H'+lambda)*(H*T(n,:)');
% % Reformulate Network
% end
fprintf('the value of Ridge Regression k')
disp(lambda)
LW=(H*H'+2*lambda*eye(size(IW,1)))\(H*Output');
%toc
%% Yalmip solver(mosek)
tic
% N=size(LW,1);
% for i = 1:S
% LWso=sdpvar(N,1);
% Obj=norm(H'*LWso-T(i,:)',2)+lambda*norm(LWso,2);
% optimize([],Obj);
% LW(:,i)=value(LWso);
% end
%Ridge Regression
toc;
fprintf('time of training ELM is %d', toc)
fprintf(' seconds \n')
fprintf('the value of Guaranteed Error of Lipschitz')
disp(norm(H'*LW-T',inf)+lambda*norm(LW,inf))
ELM1.weight{2} = LW';
X=LW'*H-T;
ELM1=ELM(ELM1.weight,ELM1.bias,ELM1.activeFcn);
ELM1.trainingError = partitions.MeanSquare(X);
end



