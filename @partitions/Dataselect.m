function [TrainDatainput,TrainDataoutput]=Dataselect(traindatainput,traindataouput,inputspace,dimension,mu,coeff,inputdimension)
%% Note:
%   Select Data according to inputspace
%
%
%inputspace=inputspace(end-1:end,:);
Num = size(inputspace,1)/dimension;
input=traindatainput';
output=traindataouput'; 
flag=1;
if (nargin<5)
 for k = 1:Num
     j=1;
     for i = 1:size(traindatainput',2)
            if(partitions.ifin(input(1:dimension,i),inputspace(dimension*(k-1)+1:dimension*k,:),dimension)==1)
                            TrainDatainput{k}(:,j)=input(:,i);
                            TrainDataoutput{k}(:,j)=output(:,i);
                            j=j+1;
                            flag=0;
            end
     end
 end
elseif(nargin==6)
for k = 1:Num
     j=1;
     for i = 1:size(traindatainput',2)
            if(partitions.ifin(coeff*(input(:,i)-mu),inputspace(dimension*(k-1)+1:dimension*k,:),dimension)==1)
                            TrainDatainput{k}(:,j)=input(:,i);
                            TrainDataoutput{k}(:,j)=output(:,i);
                            j=j+1;
                            flag=0;
            end
     end
end
elseif(nargin==7)
for k = 1:Num
     j=1;
     for i = 1:size(traindatainput',2)
            if(partitions.ifin(coeff*(input(1:inputdimension,i)-mu),inputspace(dimension*(k-1)+1:dimension*k,:),dimension)==1)
                            TrainDatainput{k}(:,j)=input(:,i);
                            TrainDataoutput{k}(:,j)=output(:,i);
                            j=j+1;
                            flag=0;
            end
     end
end
end
if(flag==1)
 TrainDatainput=[];
 TrainDataoutput=[];
end
end