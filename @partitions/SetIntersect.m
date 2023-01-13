function InterSet=SetIntersect(setinput,partitions)
% To evaluate the intersection of input set and partitions
% The dimension of setinput is (2*dimension,1)
% The dimension of one partition is (dimension,2)
dimension= size(setinput,1)/2;
Num=size(partitions,1);
InterSet=zeros(2*dimension,1);

for i = 1:Num
    for j = 1:dimension
        if partitions(j,1)>setinput(2*j-1,1)
            InterSet(2*j-1,1)= partitions(j,1);
        else
            InterSet(2*j-1,1)= setinput(2*j-1,1);
        end
        if partitions(j,2) > setinput(2*j,1)
            InterSet(2*j,1)= setinput(2*j,1);
        else
            InterSet(2*j,1)= partitions(j,2);
        end
    end
end

end