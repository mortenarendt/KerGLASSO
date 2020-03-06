function results = getEBIC(results)

% ebic
n = size(results.data{1,1},1);
p = size(results.S{1},1);
pEval = p*(p-1)/2;
n_edges = pEval  - results.n0;
eBIC05 = -n*results.loglik(2,:) + n_edges*log(n) + 0.5*4*n_edges*log(p);
eBIC00 = -n*results.loglik(2,:) + n_edges*log(n);
eBIC1 = -n*results.loglik(2,:) + n_edges*log(n) + 1*4*n_edges*log(p);

eBIC = [eBIC00; eBIC05; eBIC1];

results.ebic_gamma = [0 0.5 1];
results.eBIC = eBIC;
