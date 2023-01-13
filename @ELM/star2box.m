function outputvec= star2box(star)
Star_Num=size(star,2);
for i =1:Star_Num
    B(1,i)=getBox(star(1,i));
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
outputvec=output;
end