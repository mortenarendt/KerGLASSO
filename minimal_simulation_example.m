clc; clear all; close all; 

% Simulate a five node multivariate graph
n = [100 300];
s2n = 1;
noisetype = 1;
[A, XX, TT, G, S] = setup_simulation_data(n, s2n , [noisetype 1]);

Theta0 = zeros(size(S{1}));
fixthese = [];
max_iter = 150;
% set penalty space
PP = logspace(log10(0.001),log10(1),20)
 % DP-GLASSO
[Theta n0] = glasso_path(S{1},PP,Theta0,fixthese,'glasso_v3');

% Evaluate on test-set
loglik = getLogLikS(Theta,S{2},PP);
% Compare with known Graph
AUC = []; 
for i=1:size(Theta,1)
    th = squeeze(Theta(i,:,:));
    AUC(i) = getROC(G,th);   
end


