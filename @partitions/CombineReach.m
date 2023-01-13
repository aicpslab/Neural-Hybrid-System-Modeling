function outputbound=CombineReach(inputbound)
         for i = 1:size(inputbound,1)/2
         outputbound(2*i-1,1) = min(inputbound(2*i-1,:));
         outputbound(2*i,1) = max(inputbound(2*i,:));
         end
   
end