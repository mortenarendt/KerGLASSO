function [A, XX, TT, G, S] = setup_simulation_data(N,s2n, settings)

% Function making connected graph of 5 matrices all with similar signal to
% (white) noise and the same number of samples. Different number of variables and
% different sets of bilinear components.

% settings : Boolean (1 x 2) with 1: white (1)- or structured noise (2) and 2: remove diagonal

n = sum(N);
% n = 100; % number of individuals
k = ceil(rand(5,1)*6 + 3);
p = ceil(rand(5,1)*130+ 20);
%s2n = 10;
% make X1, X2  and X3 totally independent
% make X4 = f(X1,X2)
% make X5 = f(X2,X4)
G = eye(5);
G(4,[1 2]) = 1;
G(5,[2 4]) = 1;
% the relation is encoded in sampling sets of common "scores" and making
% arbitrary loadings and add noise

TT{1} = randn(n,k(1));
TT{2} = randn(n,k(2));
TT{3} = randn(n,k(3));
TT{4} = TT{1}*randn(k(1),k(4)) + TT{2}*randn(k(2),k(4));
TT{5} = TT{2}*randn(k(2),k(5)) + TT{4}*randn(k(4),k(5));

if settings(1)==1 % add white noise to system
    for i=1:5
        xi = TT{i}*randn(k(i),p(i));
        E = randn(size(xi));
        s2nobs = trace(xi*xi') / trace(E*E');
        E = E* sqrt(s2nobs / s2n);
        XX{i} = xi + E;
    end
elseif settings(1)==2 % add structured noise to system
    for i=1:5
        E = randn(size(TT{i}));
        s2nobs = trace(TT{i}*TT{i}') / trace(E*E');
        E = E* sqrt(s2nobs / s2n);
        xi = (TT{i} + E)*randn(k(i),p(i));
        XX{i} = xi;
    end
    
end

% split into train (N(1) and test N(2))

for i=1:size(XX,2)
    X2{1,i} = XX{i}(1:N(1),:);
    X2{2,i} = XX{i}((N(1)+1):end,:);
end

XX = X2;

% kernelize
for j=1:2
    a3 = [];
    for i=1:5
        xi = auto(XX{j,i});
        Ki = euclideandist(xi);
        mk = mean(Ki(:));
        KK{i} = exp(-Ki / mk);
        a3 = [a3 KK{i}(:)];
    end
    A{j} = a3;
end

if settings(2)==1
    for j=1:2
        iddiag = getDiagID(N(j));
        A{j} = A{j}(iddiag==0,:);
    end
end

% standardize
for j=1:2
    A{j} = A{j} - ones(size(A{j},1),1)*mean(A{j});
    A{j} = A{j}*diag(1./std(A{j}));
    S{j} = A{j}'*A{j} / size(A{j},1);
end


function iddiag = getDiagID(n)
ID = ones(n,1)*(1:n);
tID = ID';
iddiag = ID(:) - tID(:)==0;
