clear
close all
clc
% How many trajectories we want to obtain
TrajNum=4000;
NeuronNum_switch=40;
NeuronNum_single=200;
% Sampling data from working zone
u=0.2;
delta=0.5;
tau=0.01;
thetaInterval=[-4 4;-1*4 4];
duration = 150;
e=3e-5;
tol = 1;
maximum_entropy=40;
Dimension=2;
% Initial Input Set
%durationReach=20;
%ubound = [-0.08;0.02] ;
%inputbound = [2.8;3;2.5;2.6];

d = 0.1;
delta=d/2;
N = 20;
%% Robotic arm model
l1 = 10; % length of first arm
l2 = 7; % length of second arm
theta1_start = 0;
theta1_end = pi/2;
theta2_start = 0;
theta2_end = pi;
theta1 = theta1_start:d:pi/theta1_end; % all possible theta1 values
theta2 = theta2_start:d:theta2_end;
[THETA1,THETA2] = meshgrid(theta1,theta2);
X = l1 * cos(THETA1) + l2 * cos(THETA1 + THETA2); % compute x coordinates
Y = l1 * sin(THETA1) + l2 * sin(THETA1 + THETA2); % training output data
inputData = [THETA1(:),THETA2(:)]';
outputData = [X(:),Y(:)]';
% Umgebung
% theta1_u = theta1-1/2*d;
% theta1_l = theta1+1/2*d;
% theta2_u = theta2-1/2*d;
% theta2_l = theta2+1/2*d;
% [Theta1_u,Theta2_u]=meshgrid(theta1_u,theta2_u);
% [Theta1_l,Theta2_l]=meshgrid(theta1_l,theta2_l);
% Input_u = [Theta1_u(:),Theta2_u(:)]';
% Input_l = [Theta1_l(:),Theta2_l(:)]';

%
[S,~] = size(inputData);
[R,Q] = size(outputData);
%% Train ELM with robust optimization method
%1. Randomly Generate the Input Weight Matrix
    IW = rand(N,R) * 2 - 1;
%2. Randomly Generate the Bias Matrix
    B = rand(N,1);
    BiasMatrix = repmat(B,1,Q);
%3. Compute the Layer Output Matrix H
    [H_l,H_u]=computelayerOutput(Input_l',Input_u',IW,B,'sig');
    H=1/2*(H_l+H_u);
    H1=1/2*(H_u-H_l);
%4. RO SDP
   %SDP var
    LW = sdpvar(N,2);
    lamda = sdpvar(2,2);
    tao = sdpvar(2,2);
    M = H1*LW;
    F =(H*LW-outputData');
    l = norm(H*LW-outputData');
   %SDP constrain
    Cons=[tao>=0;[lamda-tao zeros(2,2) F';zeros(2,2) tao M';F M eye(size(M,1))]>=0];
   optimize(Cons,lamda);
   LW=value(LW);

X=inputData;
ELM.weight{1} = IW;
ELM.weight{2} = LW';
ELM.bias{1} = B;
ELM.bias{2} = 0;
ELM.activeFcn = {'sig','purelin'};
Y=elmpredict(X,ELM);
%% ELM
tempH = IW * inputData + BiasMatrix;
H = 1 ./ (1 + exp(-tempH));
% Calculate the Output Weight Matrix
LW1 =pinv(H')*outputData';
ELM1=ELM;
ELM1.weight{2}=LW1';
Y1=elmpredict(X,ELM1);
ERR=LW-LW1;
%% Output guaranteed distance
% Output Interval
options.tol = d*0.1;
tol=options.tol;
%Outputset of Robotic arm model 
THEta1 = theta1_start:tol:pi/theta1_end; % all possible theta1 values
THEta2 = theta2_start:tol:theta2_end;
[THETA1,THETA2] = meshgrid(THEta1,THEta2);
X = l1 * cos(THETA1) + l2 * cos(THETA1 + THETA2); % compute x coordinates
Y = l1 * sin(THETA1) + l2 * sin(THETA1 + THETA2); % training output data
outputData = [X(:),Y(:)]';
for i= 1:size(outputData,2)
      for j=1:2
           outputset{1,i}(1,1) = outputData(1,i)-l1*tol;
           outputset{1,i}(1,2) = outputData(1,i)+l1*tol;
           outputset{1,i}(2,1) = outputData(2,i)-l2*tol;
           outputset{1,i}(2,2)=  outputData(2,i)+l2*tol;
     end
end

%Outputset of ELM
inputIntvl=[0,pi/2;0,pi];
elm_ro = ffnetwork(ELM.weight,ELM.bias,ELM.activeFcn);
elm = ffnetwork(ELM1.weight,ELM1.bias,ELM1.activeFcn);
yInterval_ro = outputSet(elm_ro,inputIntvl,options);
yInterval_el = outputSet(elm,inputIntvl,options);
%1.Output mixed method for ELMs
%dis_ro=zeros(size(yInterval_ro,2),1);
%dis_el=zeros(size(yInterval_ro,2),1);
%
distance_ro=zeros(length(yInterval_ro),1);
distance_el=zeros(length(yInterval_el),1);
lossmax_ro=zeros(length(yInterval_ro),2);
lossmax_el=zeros(length(yInterval_ro),2);
for j = 1:length(yInterval_ro)
    for i= 1:2
        lossmax_ro(j,i) = max([outputset{1,j}(i,2)-yInterval_ro{1,j}(i,1),yInterval_ro{1,j}(i,2)-outputset{1,j}(i,1)]);
        lossmax_el(j,i) = max([outputset{1,j}(i,2)-yInterval_el{1,j}(i,1),yInterval_el{1,j}(i,2)-outputset{1,j}(i,1)]);
    end
end
for i=1:length(yInterval_ro)
   distance_ro(i,1)=norm(lossmax_ro(i,:)',2);
   distance_el(i,1)=norm(lossmax_el(i,:)',2);
end
e_max_reach_ro = max(distance_ro);
e_max_reach_el = max(distance_el);
%
inputData = [THETA1(:),THETA2(:)]';
y_el=elmpredict(inputData,ELM1);
y_ro=elmpredict(inputData,ELM);
%
dist_sample_ro = max(vecnorm(y_ro-outputData));
dist_sample_el = max(vecnorm(y_el-outputData));
Lipschitz_ro =dist_sample_ro+(0.25*norm(ELM.weight{1,1},2)*norm(ELM.weight{1,2},2)+sqrt(2)*(l1+l2))*tol;
Lipschitz_el =dist_sample_el+(0.25*norm(ELM1.weight{1,1},2)*norm(ELM1.weight{1,2},2)+sqrt(2)*(l1+l2))*tol;
%% Plot figures 
% figure
% hold on
% plot(X,Y,'o',X,Y1,'*',X,sin(X),'+');
% legend('ELM using RO','ELM','sin(x)')
% title('function and ELMs')
% figure
% error1=Y-Y1;
% error2=sin(X)-Y1;
% error3=sin(X)-Y;
% plot(X,error2,'*',X,error3,'o')
% legend('error between ELM and sin(x)','error between RO and sin(x)')
% title('ERRORS')

