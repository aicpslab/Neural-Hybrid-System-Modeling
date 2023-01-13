function [Y_l,Y_u]=computelayerOutput(X_l,X_u,W,bias,activFun)
[Q,R]=size(X_l);
N=size(W,1);
Bias=repmat(bias,1,Q);
P_min=X_l';
P_max=X_u';
H_l = zeros(N,Q);
H_u = zeros(N,Q);
for i = 1:Q
    for j = 1:N
        for k = 1:R
            if W(j,k) > 0
              H_l(j,i) = W(j,k)*(P_min(k,i))+H_l(j,i);
              H_u(j,i) = W(j,k)*(P_max(k,i))+H_u(j,i);
            else
              H_l(j,i)=(W(j,k)*P_max(k,i))+H_l(j,i);
              H_u(j,i)=(W(j,k)*P_min(k,i))+H_u(j,i);  
            end
        end
    end
end
% Activiation function
switch(activFun)
    case'sig'
     H_l = 1 ./ (1 + exp(-(H_l+Bias)));
     H_u = 1 ./ (1 + exp(-(H_u+Bias)));
    case'tansig'
     H_l = (1-exp(H_l+Bias))./(1+exp(H_l+Bias));   
     H_u = (1-exp(H_u+Bias))./(1+exp(H_u+Bias));
    case'purelin'
     H_l=H_l+Bias;
     H_u=H_u+Bias;
    case'ReLu'
     H_l = max(H_l+Bias,0);
     H_u = max(H_u+Bias,0);  
end
Y_l = H_l';
Y_u = H_u';
end