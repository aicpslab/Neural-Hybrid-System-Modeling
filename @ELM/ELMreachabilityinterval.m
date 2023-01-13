function output=ELMreachabilityinterval(InputBound,ELM)
% Reachability Anlynasis of ELM
% Input: 
% 1.Inputbound(2N*1): [Lowerbound_i;Upperbound_i,...] of ELM
% 2.ELM
% Output:
% 3.Output of ELM
[N,S]=size(ELM.weight{1,1});
R=size(ELM.weight{1,2},1);
inputweight=ELM.weight{1,1};
outputweight=ELM.weight{1,2};
inputbias=ELM.bias{1,1};
%% Input upper and lower bound and weight matrices 
for i = 1:size(InputBound,1)/2
    input_lower(i,1)=InputBound(2*i-1,1);
    input_upper(i,1)=InputBound(2*i,1);
end
inputweight_upper=zeros(N,S);
inputweight_lower=zeros(N,S);
outputweight_upper=zeros(R,N);
outputweight_lower=zeros(R,N);
for i = 1:N
    for j = 1:S
        if inputweight(i,j)>0
            inputweight_upper(i,j)=inputweight(i,j);
        else
            inputweight_lower(i,j)=inputweight(i,j);
        end
    end
end
for i = 1:R
    for j = 1:N
        if outputweight(i,j)>0
             outputweight_upper(i,j)=outputweight(i,j);
        else
             outputweight_lower(i,j)=outputweight(i,j);
        end
    end
end
%% Hidden layer
Hidden_lower=inputweight_upper*input_lower+inputweight_lower*input_upper+inputbias;
Hidden_upper=inputweight_lower*input_lower+inputweight_upper*input_upper+inputbias;
switch ELM.activeFcn{1,1}
    case 'sig'
        Hl = 1 ./ (1 + exp(-Hidden_lower));
        Hu = 1 ./ (1 + exp(-Hidden_upper));
    case 'sin'
        Hl = sin(Hidden_lower);
        Hu = sin(Hidden_upper);
    case 'hardlim'
        Hl = hardlim(Hidden_lower);
        Hu = hardlim(Hidden_upper);
    case 'ReLu'
        Hl = max(Hidden_lower,0);
        Hu = max(Hidden_upper,0);
end
%% Output Layer
outputlower=outputweight_upper*Hl+outputweight_lower*Hu;
outputupper=outputweight_lower*Hl+outputweight_upper*Hu;
for i =1: size(outputlower,1)
 output(2*i-1,1)=outputlower(i,1);
 output(2*i,1)=outputupper(i,1);
end
end
