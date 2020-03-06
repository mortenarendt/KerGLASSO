function dist = myeuclideandist(x)

%%
n = size(x,1); 
for i = 1:n-1
    for j = (i+1):n
        d = x(i,:) - x(j,:); 
        sq = sqrt(d*d'); 
        dist(i,j) = sq;
        dist(j,i) = sq;
    end
end


