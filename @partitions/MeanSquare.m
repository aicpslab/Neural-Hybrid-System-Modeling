function Y=MeanSquare(X)
[N,Q]=size(X);
y=0;
for i = 1:Q
    for j = 1:N
    y=y+X(j,i)^2;
    end
end
Y=1/Q*(sqrt(y));
end