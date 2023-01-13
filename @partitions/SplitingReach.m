function [splitingset,t] = SplitingReach(intervals,ELM,setinput,Reachmethod,ubound)
ELMbound=setinput;
inputspace1=intervals;
counter=1;
dimension = size(setinput,1)/2;

tic


            for k = 1:size(inputspace1,2)
               for z = 1:size(inputspace1{k},1)/dimension
                   inputbound=partitions.SetIntersect(setinput,inputspace1{k}((dimension*(z-1)+1:dimension*z),:));
%                 if inputspace1{1,k}(2*z-1,1)>ELMbound(1,1)
%                     inputbound(1,1)= inputspace1{1,k}(2*z-1,1);
%                 else
%                     inputbound(1,1)= ELMbound(1,1);
%                 end
%                 if inputspace1{1,k}(2*z,1) > ELMbound(3,1)
%                     inputbound(3,1)= inputspace1{1,k}(2*z,1);
%                 else
%                     inputbound(3,1)= ELMbound(3,1);
%                 end
%                 if inputspace1{1,k}(2*z-1,2) < ELMbound(2,1)
%                     inputbound(2,1)= inputspace1{1,k}(2*z-1,2);
%                 else
%                     inputbound(2,1)= ELMbound(2,1);
%                 end
%                 if inputspace1{1,k}(2*z,2) < ELMbound(4,1)
%                     inputbound(4,1)= inputspace1{1,k}(2*z,2);
%                 else
%                     inputbound(4,1)= ELMbound(4,1);
%                 end
                if(nargin<=4)
                    ubound=[];
                end
                if (inputbound(1,1)<inputbound(2,1)&&inputbound(3,1)<inputbound(4,1))
                    InputBound(:,counter) = [inputbound(:,1);ubound];
                switch(Reachmethod)
                     case 'nnv'
                       outputbound(1:4,counter) = ELMreachabilitynnv(InputBound(:,counter),ELM(k));
                     case 'interval'
                       outputbound(1:4,counter) = ELMreachabilityinterval(InputBound(:,counter),ELM(k));
                end
                    counter=counter+1;
                end
               end
           end
splitingset=outputbound;
t=toc;
end