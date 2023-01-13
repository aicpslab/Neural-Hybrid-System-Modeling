clear
close all
clc
%load(['DataSet/inverse_identification_without_raw_data.mat']) %loading the data set

% How many trajectories we want to obtain
NeuronNum_switch=40;
NeuronNum_single=400;
% Sampling data from working zone

TF='sig';
%duration = size(y_test,2);
duration =3636;

e=4e-4;


tol = 1.5;
maximum_entropy=40;

% Initial Input Set
%durationReach=20;

%% Load data set
load(['DataSet/forward_identification_without_raw_data.mat'])
xs = zeros(size(u_train,2)-1,18);
t = zeros(size(u_train,2)-1,6);
for i = 1:size(u_train,2)-2
   xs(i,:)=[y_train(:,i)',y_train(:,i+1)',u_train(:,i+1)'];
   t(i,:)=[y_train(:,i+2)'];  
end
for i = 1:size(y_test,1)
    subplot(size(y_test,1),1,i)
    plot(time_test,y_test(i,:))
    xlabel('time(s)')
    ylabel('position')
    hold on
end




delta_y = zeros(6,size(time_train,2)-1);

for i = 1:size(time_train,2)-1
   delta_y(:,i)=abs(y_train(:,i+1)-y_train(:,i));
end
MaxDelta_y=zeros(size(y_train,1),1);
for i =1:6
 MaxDelta_y(i,1)=max(delta_y(i,:));
end
%xs=u_train';
%t=y_train';

[Pn_train,inputps] = mapminmax(xs');
%Pn_test = mapminmax('apply',P_test,inputps);
[Tn_train,outputps] = mapminmax(t');
%% Princeple Component Analysis
% [coeff,scoreeeTrain,~,~,explained,mu] = pca(xs(:,1:12));
% idx=find(cumsum(explained)>90,1);      
% utest = (xs(:,1:12)-mu)*coeff(:,1:idx);
% bounderies = zeros(idx,2);
[coeff,scoreeeTrain,~,~,explained,mu] = pca(Pn_train(1:12,:)');
idx=find(cumsum(explained)>90,1);      
utest = (Pn_train(1:12,:)'-mu)*coeff(:,1:idx);
bounderies = zeros(idx,2);

for i= 1:idx
    bounderies(i,1)=min(utest(:,i))-2;
    bounderies(i,2)=max(utest(:,i))+2;
end
lowerbound=bounderies(:,1);
upperbound= bounderies(:,2);
%init_interval{1}=[lowerbound',upperbound'];
init_interval{1}=bounderies;
Dimension=size(bounderies,1);

%% Data-driven Partitioning
P=partitions(init_interval,Pn_train(1:12,:)',Tn_train');
intervals=ME(P,tol,maximum_entropy,Dimension,mu',coeff(:,1:idx)');
%intervals=ME(P,tol,maximum_entropy,Dimension,mu',coeff(:,1:idx)');
P.intervals=intervals;
P.input=Pn_train';

%% Initialized ELM and Merge Partitions
ELMs1=ELM.GenerateELM(size(xs,2),NeuronNum_switch,TF,size(t,2));
[P1,ELMs]=MergePatitions(P,ELMs1,e,mu',coeff(:,1:idx)',12);
size(P.intervals,2)
size(P1.intervals,2)
mse_switch = 0;
min_mse_switch=1;
for i = 1:size(ELMs,2)
    if (mse_switch<ELMs(i).trainingError)
            mse_switch = ELMs(i).trainingError;
            k=i;
    end
    if (min_mse_switch>ELMs(i).trainingError)
            min_mse_switch = ELMs(i).trainingError;
    end
end

%  for i = 1:size(ELMs,2)
%     ELMs(i).trainingError
%      i         
%  end

% Plot intervals
% figure
% partitions.intervalplot(P.intervals,'empty','blue')
% partitions.intervalplot(P1.intervals,'full','red')
% title('Invariantspace using Bisection method')

%% Train a Complex Neural Network Model as referance
ELMs1=ELM.GenerateELM(size(xs,2),NeuronNum_single,TF,size(t,2));
%ELMs1=trainELMLipridge(ELMs1,xs',t');
ELMs1=trainELM(ELMs1,Pn_train,Tn_train);
%ELMs1=trainELMLipridge(ELMs1,Pn_train,Tn_train);

%% Verfify whether can approximate the model well
output_single=zeros(6,size(y_test,2)-2);
output_switch=output_single;
segmentIndex=P1.intervals;
inputspace1=P1.intervals;
%   xs(i,:)=[y_train(:,i)',y_train(:,i+1)',u_train(:,i+1)'];
%   t(i,:)=[y_train(:,i+2)'];  

for i =1:size(time_test,2)-2
    IN=mapminmax('apply',[y_test(:,i);y_test(:,i+1);u_test(:,i+1)],inputps);
     output_single(:,i) = ELMpredict(ELMs1,IN);
  
   for k = 1:size(segmentIndex,2)
              if(partitions.ifin(coeff(:,1:idx)'*(IN(1:12,:)-mu'),segmentIndex{k},Dimension)==1)
                       output_switch(:,i)= ELMpredict(ELMs(k),IN);

              end
    end
end

%BPs=newff(xs',t',[18 200 200 6],{'purelin' 'tansig' 'tansig' 'purelin'});
%BPs=train(BPs,xs',t');
%% Analysis from the stastical point of view
 % 1.Plot Figures
figure
ON= mapminmax('reverse',output_single(:,1:end),outputps);
 for i = 1:size(y_test,1)-1
    subplot(size(y_test,1),1,i)
    plot(time_test(3:end),y_test(i,3:end))
    hold on 
   
    plot(time_test(3:end),ON(i,1:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
 end
 title('single neural network modeling')

ON_S= mapminmax('reverse',output_switch(:,1:end),outputps);
figure
for i = 1:size(y_test,1)-1
    subplot(size(y_test,1),1,i)
    plot(time_test(3:end),y_test(i,3:end))
    hold on 
    plot(time_test(3:end),ON_S(i,1:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
end
title('switching neural network modeling')

%% Simulate the Respones
Traj=zeros(size(y_test,1),duration);
Traj(:,1:2)=[y_test(:,1),y_test(:,2)];
Traj_single=Traj;
Traj_switch=Traj;
%inputN(:,1:3)= Traj

for i =3:duration
    InputN=mapminmax('apply',[Traj_single(:,i-2);Traj_single(:,i-1);u_test(:,i-1)],inputps);
    Traj_singleN=ELMpredict(ELMs1,InputN);
    Traj_single(:,i)=mapminmax('reverse',Traj_singleN,outputps);
    output_delta=Traj_single(:,i)-Traj_single(:,i-1);
   
     for j= 1:size(output_single,1)  
      if(abs(output_delta(j,1))>MaxDelta_y(j,1))
              output_delta(j,1)=MaxDelta_y(j,1)*sign(output_delta(j,1));
      end
    end  
    
    Traj_single(:,i)=Traj_single(:,i-1)+output_delta;
    
    InputNs=mapminmax('apply',[Traj_switch(:,i-2);Traj_switch(:,i-1);u_test(:,i-1)],inputps);
    
    for k = 1:size(segmentIndex,2)
              if(partitions.ifin(coeff(:,1:idx)'*(InputNs(1:12,:)-mu'),segmentIndex{k},Dimension)==1)
                       Traj_switchN= ELMpredict(ELMs(k),InputNs);
                       Traj_switch(:,i)=mapminmax('reverse',Traj_switchN,outputps);
              end
    end
    output_delta=Traj_switch(:,i)-Traj_switch(:,i-1);
    for j= 1:size(output_single,1)  
      if(abs(output_delta(j,1))>MaxDelta_y(j,1))
              output_delta(j,1)=MaxDelta_y(j,1)*sign(output_delta(j,1));
      end
    end  
     
    Traj_switch(:,i)=Traj_switch(:,i-1)+output_delta;
end
% 1.Plot Figures
figure
 for i = 1:size(y_test,1)-1
    subplot(size(y_test,1),1,i)
    plot(time_test(3:duration),y_test(i,3:duration))
    hold on 
    plot(time_test(3:duration),Traj_single(i,3:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
 end
 title('single neural network modeling')
figure
for i = 1:size(y_test,1)-1
    subplot(size(y_test,1),1,i)
    plot(time_test(3:duration),y_test(i,3:duration))
    hold on 
    plot(time_test(3:duration),Traj_switch(i,3:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
end
title('switching neural network modeling')

% 3.Linear model
% define ssest options
%  opt_ssest = ssestOptions ();
%  opt_ssest.Display = 'on';
%  opt_ssest.Focus = 'simulation';
% 
%  % zero initialize the input
%  utrain = u_train - u_train (:,1);
%  utest = u_test-u_test(:,1);
%  % define data
%  dt = 0.1;
%  id_data = iddata(y_train',utrain',dt);
%  % identify model
%  n_states = 12;
%  ss_model = ssest (id_data,n_states,opt_ssest );
%  compare(id_data,ss_model)
%  test_data=iddata(y_test',utest',dt);
%  figure
%  [y,fit,ic]=compare(test_data,ss_model)
%  compare(test_data,ss_model)

% % 4. Compare from the stastical point of view
% devia=zeros(size(Traj_single,1),1);
%  %devia=std(y_test');
% for i = 1:size(devia,1)
%   devia(i,1)=std(y_test(i,:));
% end
% NRMSE_single=zeros(size(devia,1),1);
% NRMSE_switch=zeros(size(devia,1),1);
% ERR_single=zeros(size(devia,1),1);
% ERR_switch=zeros(size(devia,1),1);
% for i = 1:size(devia,1)
%     for j = 1:size(Traj_single,2)
%         ERR_single(i,1)=ERR_single(i,1)+(Traj_single(i,j)-mean(y_test(i,:)))^2;
%         ERR_switch(i,1)=ERR_switch(i,1)+(Traj_switch(i,j)-mean(y_test(i,:)))^2;
%     end
% end
% for i = 1:size(devia,1)
%   NRMSE_single(i,1)=sqrt(1/size(Traj_single,2)*1/(devia(i,1))^2*ERR_single(i,1));
%   NRMSE_switch(i,1)=sqrt(1/size(Traj_single,2)*1/(devia(i,1))^2*ERR_switch(i,1));
% end

%% Reachable set computation for our model 
% Ini_Input_u=[y_test(:,1);y_test(:,2)];
% Ini_Input_l=[y_test(:,1);y_test(:,2)];
% 
% Ini_SetInput=zeros(2*size(Ini_Input_l,1),1);
% 
% for i = 1:size(y_test,1)
%     Ini_Input_u(size(y_test,1)+i)=Ini_Input_l(size(y_test,1)+i)+delta;
% end
% 
% for i = 1:size(Ini_Input_l,1)
%     Ini_SetInput(2*(i-1)+1,1)= Ini_Input_l(i);
%     Ini_SetInput(2*i,1)=Ini_Input_u(i);
% end
% % Input torque: u_test(:,2)
% % Bound of Input torque
% ubound = zeros(size(u_test,1)*2,size(u_test,2));
% for i = 1:size(u_test,2)
%     for j = 1:size(u_test,1)
%         ubound(2*(j-1)+1,i)=u_test(j,i);
%         ubound(2*(j),i)=u_test(j,i);
%     end
% end
% %1. Single ELM model
% fprintf('Reachable set estimation for Single ELM using NNV.\n')
% 
% singleELMtime=zeros(durationReach,1);
% for i =1:1 % How many trajectories what to be examined
%     SELMbound{i}(:,1)= Ini_SetInput(1:2*size(y_test,1),:);
%     SELMbound{i}(:,2)= Ini_SetInput(2*size(y_test,1)+1:end,:);
%     for j = 3:durationReach
%         fprintf('step'); 
%        disp(j); 
%         tic
%         inputboundELM = ELMreachabilitynnv([SELMbound{i}(:,j-2);SELMbound{i}(:,j-1);ubound(:,j-1)],ELMs1);    
%         SELMbound{i}(:,j)=[inputboundELM];
%      toc
%      singleELMtime(j-1,1)=toc;
%     end
% end
% 
% switchELMtime=zeros(durationReach,1);
% %2. Switching ELM model  
% fprintf('Reachable set estimation for Single ELM using NNV.\n')
% for i=2:durationReach
% tic
% 
% switchELMtime(j-1,1)=toc;
% end



Traj_train=zeros(size(y_train,1),size(y_train,2));
Traj_train(:,1:2)=[y_train(:,1),y_train(:,2)];
Traj_singlet=Traj_train;
Traj_switcht=Traj_train;
%inputN(:,1:3)= Traj

for i =3:size(y_train,2)
    InputN=mapminmax('apply',[Traj_singlet(:,i-2);Traj_singlet(:,i-1);u_train(:,i-1)],inputps);
    Traj_singleN=ELMpredict(ELMs1,InputN);
    Traj_singlet(:,i)=mapminmax('reverse',Traj_singleN,outputps);
    output_delta=Traj_singlet(:,i)-Traj_singlet(:,i-1);
   
     for j= 1:size(output_single,1)  
      if(abs(output_delta(j,1))>MaxDelta_y(j,1))
              output_delta(j,1)=MaxDelta_y(j,1)*sign(output_delta(j,1));
      end
    end  
    
    Traj_singlet(:,i)=Traj_singlet(:,i-1)+output_delta;
    
    InputNs=mapminmax('apply',[Traj_switcht(:,i-2);Traj_switcht(:,i-1);u_train(:,i-1)],inputps);
    
    for k = 1:size(segmentIndex,2)
              if(partitions.ifin(coeff(:,1:idx)'*(InputNs(1:12,:)-mu'),segmentIndex{k},Dimension)==1)
                       Traj_switchN= ELMpredict(ELMs(k),InputNs);
                       Traj_switcht(:,i)=mapminmax('reverse',Traj_switchN,outputps);
              end
    end
    output_delta=Traj_switcht(:,i)-Traj_switcht(:,i-1);
    for j= 1:size(output_single,1)  
      if(abs(output_delta(j,1))>MaxDelta_y(j,1))
              output_delta(j,1)=MaxDelta_y(j,1)*sign(output_delta(j,1));
      end
    end  
     
    Traj_switcht(:,i)=Traj_switcht(:,i-1)+output_delta;
end
% 1.Plot Figures
figure
 for i = 1:size(y_train,1)-1
    subplot(size(y_train,1),1,i)
    plot(time_train(3:end),y_train(i,3:end))
    hold on 
    plot(time_train(3:end),Traj_singlet(i,3:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
 end
 title('single neural network modeling')
figure
for i = 1:size(y_train,1)-1
    subplot(size(y_train,1),1,i)
    plot(time_train(3:end),y_train(i,3:end))
    hold on 
    plot(time_train(3:end),Traj_switcht(i,3:end))
    xlabel('time(s)')
    ylabel('position')
    hold on
end
title('switching neural network modeling')

devia=zeros(size(Traj_singlet,1),1);
 %devia=std(y_test');
% for i = 1:size(devia,1)
%   deviat(i,1)=std(y_train(i,:));
% end
% NRMSE_singlet=zeros(size(deviat,1),1);
% NRMSE_switcht=zeros(size(deviat,1),1);
% ERR_single=zeros(size(deviat,1),1);
% ERR_switch=zeros(size(deviat,1),1);
% for i = 1:size(deviat,1)
%     for j = 1:size(Traj_singlet,2)
%         ERR_single(i,1)=ERR_single(i,1)+(Traj_singlet(i,j)-mean(y_train(i,:)))^2;
%         ERR_switch(i,1)=ERR_switch(i,1)+(Traj_switch(i,j)-mean(y_train(i,:)))^2;
%     end
% end
% for i = 1:size(devia,1)
%   NRMSE_singlet(i,1)=sqrt(1/size(Traj_singlet,2)*1/(deviat(i,1))^2*ERR_single(i,1));
%   NRMSE_switcht(i,1)=sqrt(1/size(Traj_singlet,2)*1/(deviat(i,1))^2*ERR_switch(i,1));
% end

NRMSE_singlet=zeros(size(devia,1),1);
NRMSE_switcht=zeros(size(devia,1),1);
ERR_single=zeros(size(devia,1),1);
ERR_switch=zeros(size(devia,1),1);
x2=zeros(size(devia,1),1);
for i = 1:size(devia,1)
    for j = 1:size(Traj_singlet,2)
        ERR_single(i,1)=ERR_single(i,1)+(Traj_singlet(i,j)-y_train(i,j))^2;
        ERR_switch(i,1)=ERR_switch(i,1)+(Traj_switcht(i,j)-y_train(i,j))^2;
        x2(i,1)=x2(i,1)+y_train(i,j)^2;  
    end
end
for i = 1:size(devia,1)
  NRMSE_singlet(i,1)=sqrt(1/x2(i,1)*ERR_single(i,1));
  NRMSE_switcht(i,1)=sqrt(1/x2(i,1)*ERR_switch(i,1));
end


NRMSE_single=zeros(size(devia,1),1);
NRMSE_switch=zeros(size(devia,1),1);
ERR_single=zeros(size(devia,1),1);
ERR_switch=zeros(size(devia,1),1);
R2single=zeros(size(devia,1),1);
R2switch=zeros(size(devia,1),1);
x2=zeros(size(devia,1),1);
for i = 1:size(devia,1)
    for j = 1:size(Traj_single,2)
        ERR_single(i,1)=ERR_single(i,1)+(Traj_single(i,j)-y_test(i,j))^2;
        ERR_switch(i,1)=ERR_switch(i,1)+(Traj_switch(i,j)-y_test(i,j))^2;
        R2single(i,1)=R2single(i,1)+(y_test(i,j)-mean(Traj_single(i,:)))^2;
        R2switch(i,1)=R2switch(i,1)+(y_test(i,j)-mean(Traj_switch(i,:)))^2;
        x2(i,1)=x2(i,1)+y_test(i,j)^2;  
    end
end
for i = 1:size(devia,1)
  NRMSE_single(i,1)=sqrt(1/x2(i,1)*ERR_single(i,1));
  NRMSE_switch(i,1)=sqrt(1/x2(i,1)*ERR_switch(i,1));
end

% R2 performance
for i = 1:6
    R2single(i,1)= (1-ERR_single(i,1)/R2single(i,1));
    R2switch(i,1)= (1-ERR_switch(i,1)/R2switch(i,1));
end

disp('NRMSE_single')
mean(NRMSE_single)
disp('NRMSE_switch')
mean(NRMSE_switch)
disp('NRMSE_singlet')
mean(NRMSE_singlet)
disp('NRMSE_switcht')
mean(NRMSE_switcht)

N=size(Traj_switcht,2);
for i = 1:6
T_sim = Traj_switcht(i,:);
T_test = y_train(i,:);
R2(i) = (N * sum(T_sim .* T_test) - sum(T_sim) * sum(T_test))^2 / ((N * sum((T_sim).^2) - (sum(T_sim))^2) * (N * sum((T_test).^2) - (sum(T_test))^2));
end
N=size(Traj_switch,2);
for i = 1:6
T_sim = Traj_switch(i,:);
T_test = y_test(i,:);
R2(i) = (N * sum(T_sim .* T_test) - sum(T_sim) * sum(T_test))^2 / ((N * sum((T_sim).^2) - (sum(T_sim))^2) * (N * sum((T_test).^2) - (sum(T_test))^2));
end