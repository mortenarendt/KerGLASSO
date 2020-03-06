function AUC = getROC(G,th)
%%
tol = 1e-5;
ID = tril(ones(size(G)),-1);
% remove diagonal and upper part
t = G(ID==1);
y = abs(th(ID==1));
y = (y>tol)+0;
%
[falseAlarmRate, detectionRate, TT,  AUC] = perfcurve(t,y,1);
