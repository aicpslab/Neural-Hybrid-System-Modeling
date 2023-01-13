function flag=ifin(data,space,dimension)

exit_flag=0;
for j = 1:size(space,1)/dimension
flag=1;    
    for i= 1:dimension
       %((input(1,i)>inputspace(2*(k-1)+1,1))&&(input(1,i)<inputspace(2*(k-1)+1,2)))
       if (flag==1)  
           if (data(i,1)<space((j-1)*dimension+i,1))||(data(i,1)>space((j-1)*dimension+i,2)) 
               flag=0;
           end
       end
    end
    if(flag==1)
        exit_flag=1;
    end
end
flag=exit_flag;
end