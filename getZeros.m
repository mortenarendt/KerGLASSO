function n0 = getZeros(S)
tol = 1e-5; 
p = size(S,1); 
S = abs(S - triu(S));
n0 = sum(S(:)<tol) - p*(p+1)/2;