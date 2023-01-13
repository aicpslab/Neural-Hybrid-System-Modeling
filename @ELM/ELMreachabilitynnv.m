function output=ELMreachabilitynnv(InputBound,ELM)
% Reachability Anlynasis of ELM
% Input: 
% 1.Inputbound(2N*1): [Lowerbound_i;Upperbound_i,...] of ELM
% 2.ELM
% Output:
% 3.Output of ELM
%% Formulating nnv struct
[N,S]=size(ELM.weight{1,1});
R=size(ELM.weight{1,2},1);
inputweight=ELM.weight{1,1};
outputweight=ELM.weight{1,2};
inputbias=ELM.bias{1,1};
outputbias=ELM.bias{1,2};
L1=LayerS(inputweight,inputbias,'poslin');
L2=LayerS(outputweight,outputbias,'purelin');
ELMnnv=FFNNS([L1 L2]);
%% Formulating input interval to Star set
dimension=size(InputBound,1)/2;
lb = zeros(dimension,1);
ub = zeros(dimension,1);
for i = 1:dimension
    lb(i,1)=InputBound(2*i-1,1);
    ub(i,1)=InputBound(2*i,1);
end
%lb = [InputBound(1,1);InputBound(3,1);InputBound(5,1)];
%ub = [InputBound(2,1);InputBound(4,1);InputBound(6,1)];
%Star set input
I = Star(lb,ub);
Bound=zeros(R,2);
cores=4;
% Compute the reachable sets with exact-star set method
%[R1 , t1]= ELMnnv.reach(I,'exact-star',cores);
[R1 , t1]= ELMnnv.reach(I,'exact-star',cores);
Star_Num=size(R1,2);
for i =1:Star_Num
    B(1,i)=getBox(R1(1,i));
    if i ==1
    Bound(:,1)=B(1,i).lb;
    Bound(:,2)=B(1,i).ub;
    else
        for j = 1:R
             if(B(1,i).lb(j,1)<Bound(j,1))
               Bound(j,1) = B(1,i).lb(j,1);
             end
             if(B(1,i).ub(j,1)>Bound(j,2))
               Bound(j,2) = B(1,i).ub(j,1);  
             end
        end    
    end   
end
output = zeros(2*R,1);
for i = 1:R
    output(2*(i-1)+1)=Bound(i,1);
    output(2*i)=Bound(i,2);
end 

%output=[Bound(1,:)';Bound(2,:)'];
%% Note that its a star set and should be converted into interval set 
end
