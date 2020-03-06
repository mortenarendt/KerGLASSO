function [aX mX sX] = myauto(X)

% calculate mean
mX = mean(X);
%%
% split data and do scaling
[n v] = size(X);

k = 10000;
c = ceil(v/k);
for i=1:c;
    st = (i-1)*k+1;
    slut = k*i;
    if i==c;
        slut = v;
    end
    id{i} = st:slut;
end

aX = zeros(n,v);
sX = zeros(1,v);
for i=1:c
    x = X(:,id{i});
    sx = std(x);
    mx = mean(x); 
    sX(id{i}) = sx;
    sx = ones(n,1)*sx;
    ax = x.*sx;
    
    ax = x - ones(n,1)*mx;
    
    aX(:,id{i}) = ax;
end






