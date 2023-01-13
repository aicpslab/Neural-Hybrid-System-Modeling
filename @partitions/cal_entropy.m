function  delta_entropy= cal_entropy(inputspace1,inputspace2,P,T,dimension)
%Caculate the interval of the outputset
%-Input £º
%
% P    - Input of the Training set
% P(R*Q) T(S*Q)
% inputspace(R,2) % Pertubation interval

% 1.Caculate the number of samples in inputspace1 and inputspace2

N = size(P,1);
[Xin_1,~]= partitions.Dataselect(P,T,inputspace1,dimension); 
[Xin_2,~]= partitions.Dataselect(P,T,inputspace2,dimension);
N1=size(Xin_1{1},2);
N2=size(Xin_2{1},2);

% 2.Caculate the viariation in the entropy
delta_entropy = (N2*log2(1+N1/N2)+N1*log2(1+N2/N1))/N;
end