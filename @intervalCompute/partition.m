function I_array = partition(I,tol)
% partition an interval
%
% Syntax:
%    I_array = partition(I,tol)
%
% Inputs:
%    I - interval matrix
%    tol - tolerance for terminating partitioning
%
% Outputs:
%    I_array - an array of interval matrices after partitioning

% Author:       Weiming Xiang
% Written:      02/25/2019
% Last update:  02/25/2019

%------------- BEGIN CODE --------------
% M{1} = I;
% numM = 1;
% numI_array = 0;
% while numM > 0
%     tempM = M{1};
%     M(1) = [];
%     if intervalCompute.terminate(tempM,tol) == 0
%         [I1,I2] = intervalCompute.bisect(tempM);
%         M{numM} = I1;
%         numM = numM+1;
%         M{numM} = I2;
%     elseif intervalCompute.terminate(tempM,tol) == 1
%         numI_array = numI_array+1;
%         I_array{numI_array} = tempM;
%         numM = numM-1;
%     end
% end

[dim,~] = size(I);
for i = 1:1:dim
    if I(i,2) - I(i,1) <= tol        
        x_dim{i} = [I(i,1),I(i,2)];
    else
        x_dim{i} = I(i,1):tol:I(i,2);
        if x_dim{i}(end) ~= I(i,2)
           x_dim{i} = [x_dim{i}, I(i,2)];
        end
    end
    for j = 1:1:length(x_dim{i})-1
        I_dim{i}{j} = [x_dim{i}(j),x_dim{i}(j+1)];
    end
end

t0 = ['k=0; '];
t1 = ['for i = 1:1:length(I_dim{1}) I1 = I_dim{1}{i}; ',];
t_end = ['end;'];
for j = 2:1:dim
    t2 = ['for i = 1:1:length(I_dim{', num2str(j),'})',' I', num2str(j),' = I_dim{',num2str(j),'}{i}; '];
    t1 = [t1 t2];
end
t3 = ['k = k+1; I_array{k} = ['];
for j = 1:1:dim
    t3 = [t3 'I' num2str(j)  ';'];
end
t3 = [t3 ']; '];
for j = 2:1:dim
    t_end = [t_end ' end;'];
end
t = [t0 t1 t3 t_end];
eval(t)

%% use built-in function ndgrid


%%
% [dim,~] = size(I);
% for i = 1:1:dim
%     if I(i,2) - I(i,1) < tol
%         x_dim{i} = [I(i,1),I(i,2)];
%     else
%         x_dim{i} = I(i,1):tol:I(i,2);
%     end
%     for j = 1:1:length(x_dim{i})-1
%         I_dim{i}{j} = [x_dim{i}(j),x_dim{i}(j+1)];
%     end
% end
% t1 = ['k=0; '];
% t2 = ['i(1) = 0; '];
% for j = 2:1:dim
%     t2 = [t2 'i(',num2str(j),') = 0; '];
% end
% t3 = ['while i(1) <= length(I_dim{1}) i(1) = i(1) + 1; I1 = I_dim{1}{i(1)}; '];
% for j = 2:1:dim
%     t3 = [t3 'while i(',num2str(j), ') <= length(I_dim{',num2str(j), '}) '];
%     t3 = [t3 'i(',num2str(j),') = i(',num2str(j),') + 1; '];
%     t3 = [t3 'I',num2str(j),' = I_dim{',num2str(j),'}{i(',num2str(j),')}; '];
% end
% t4 = ['k = k+1; I_array{k} = ['];
% for j = 1:1:dim
%      t4 = [t4 'I' num2str(j)  ';'];
% end
% t4 = [t4 ']; '];
% 
% for j = 1:1:dim
%     q = dim+1-j;
%     t4 = [t4 'if i(',num2str(q),') == length(I_dim{',num2str(q),'}) i(',num2str(q),') = 0; break;  end; end; '];
% end
% 
% t = [t1 t2 t3 t4];
% eval(t);

end


%------------- END OF CODE --------------