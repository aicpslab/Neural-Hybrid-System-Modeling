function partitions=ME(obj,length,epsilon,dimension,mu,coeff)
% partitions through Maximum Entropy method with threshold epsilon
% Syntax: partitions= obj.ME(obj,length,epsilon)
% Input: 
% obj partitions and input and target output matrix
% Note: The ME method is an abstract method which for each partition input can be
% make of multiple intervals
k=1;
if nargin < 4
    error('MEerror:Arguments','Not enough input arguments.');
end
 for i = 1:1:size(obj.intervals,1)
       for j =0:dimension:(size(obj.intervals{i},1)-dimension)
         M_X(dimension*k-dimension+1:dimension*k,:) = obj.intervals{i}(j+1:j+dimension,:);
         k=k+1;
       end
end
xs=obj.input;
N=size(xs,1);
%t=obj.output;
  tempxs=xs;
% tempt=t;
fprintf('Starting Maximum Entropy Partitioning with epsilon=%d',epsilon);
fprintf('\n');
%fprintf(/n)
k=1; 
flag=1;
maximum_entropy = epsilon;

if nargin < 5
    mu=0;
    coeff=0;
while (flag==1)
    fprintf('Number of Partitions %d',size(M_X,1)/dimension);
    fprintf('\n');
    tic
        Xtemp=M_X(dimension*k-dimension+1:dimension*k,:);
        max(abs(M_X(:,1)-M_X(:,2)));
        [value,position] = max(abs(Xtemp(:,1)-Xtemp(:,2)));             
        if value>length
                X_1 = Xtemp(position,1)+value/2;
                X_2 = Xtemp(position,2)-value/2;        
                Xtemp1 = Xtemp;
                Xtemp1(position,1) = X_1;
                Xtemp2 = Xtemp;
                Xtemp2(position,2) = X_2;
           % tic
              [input1,~]=obj.Dataselect(obj.input,obj.output,Xtemp1,dimension);
              [input2,~]=obj.Dataselect(obj.input,obj.output,Xtemp2,dimension);
           %toc
              
                if(~isempty(input1))&&(~isempty(input2))
                      N1=size(input1{1},2);
                      N2=size(input2{1},2);
                      delta_entropy = (N2*log2(1+N1/N2)+N1*log2(1+N2/N1))/N;
                    %if (obj.cal_entropy(Xtemp1,Xtemp2,obj.input,obj.output,dimension)*size(tempxs,1)>maximum_entropy)
                    if (delta_entropy*size(tempxs,1)>maximum_entropy)
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp1;
                        M_X(size(M_X,1)+1:size(M_X,1)+dimension,:) = Xtemp2;
                      %  xs=[input1{1,1},input2{1,1}];
                    else
                        k=k+1;
                  %   xs=tempxs;
                  %   t=tempt;
                    end
                elseif(~isempty(input1))
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp1;
                   %     xs=input1{1,1};
                elseif(~isempty(input2)) 
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp2;
                     %   xs=input2{1,1};
                    %    
                end
        else
            k=k+1;
        end

    if(k==size(M_X,1)/dimension+1)
       flag=0;
    end
    toc    
end
   
else

while (flag==1)
    fprintf('Number of Partitions %d',size(M_X,1)/dimension);
    fprintf('\n');
    tic
        Xtemp=M_X(dimension*k-dimension+1:dimension*k,:);
        max(abs(M_X(:,1)-M_X(:,2)));
        [value,position] = max(abs(Xtemp(:,1)-Xtemp(:,2)));             
        if value>length
                X_1 = Xtemp(position,1)+value/2;
                X_2 = Xtemp(position,2)-value/2;        
                Xtemp1 = Xtemp;
                Xtemp1(position,1) = X_1;
                Xtemp2 = Xtemp;
                Xtemp2(position,2) = X_2;
           % tic
              [input1,~]=obj.Dataselect(obj.input,obj.output,Xtemp1,dimension,mu,coeff);
              [input2,~]=obj.Dataselect(obj.input,obj.output,Xtemp2,dimension,mu,coeff);
           %toc
              
                if(~isempty(input1))&&(~isempty(input2))
                      N1=size(input1{1},2);
                      N2=size(input2{1},2);
                      delta_entropy = (N2*log2(1+N1/N2)+N1*log2(1+N2/N1))/N;
                    %if (obj.cal_entropy(Xtemp1,Xtemp2,obj.input,obj.output,dimension)*size(tempxs,1)>maximum_entropy)
                    if (delta_entropy*size(tempxs,1)>maximum_entropy)
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp1;
                        M_X(size(M_X,1)+1:size(M_X,1)+dimension,:) = Xtemp2;
                      %  xs=[input1{1,1},input2{1,1}];
                    else
                        k=k+1;
                  %   xs=tempxs;
                  %   t=tempt;
                    end
                elseif(~isempty(input1))
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp1;
                   %     xs=input1{1,1};
                elseif(~isempty(input2)) 
                        M_X(dimension*k-dimension+1:dimension*k,:)= Xtemp2;
                     %   xs=input2{1,1};
                    %    
                end
        else
            k=k+1;
        end

    if(k==size(M_X,1)/dimension+1)
       flag=0;
    end
    toc    
end

end

 for i=1:1:size(M_X,1)/dimension
        partitions{i}=M_X(dimension*i-dimension+1:dimension*i,:);
 end


end