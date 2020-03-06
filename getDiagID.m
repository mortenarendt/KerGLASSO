function iddiag = getDiagID(n)
ID = ones(n,1)*(1:n);
tID = ID';
iddiag = ID(:) - tID(:)==0;
end
