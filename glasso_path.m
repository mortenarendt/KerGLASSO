function [Theta, n0,round] = glasso_path(S_train,L,Theta0,fixthese,method)

if nargin==2
    Theta0 = [];
    fixthese = [];
end

Theta = nan([length(L), size(S_train)]);

for i=1:length(L)
    % train model
    switch method
        %%
        %clc
        %ll = 0.05
        %Theta_v2    = glasso_v2(S_train,ll,Theta0,fixthese);
        %Theta_v3    = glasso_v3(S_train,ll,Theta0,fixthese);
        %Theta_org   = glasso(   S_train,ll,Theta0,fixthese);
        %Theta_v3(abs(Theta_v3)<1e-5) = 0;
        %corr([Theta_v2(:),Theta_v3(:),Theta_org(:)])
        %%
        case 'glasso_v3'
            [Theta(i,:,:) ,~,round(i)] = glasso_v3(   S_train,L(i),Theta0,fixthese);
        case 'glasso_v2'
            Theta(i,:,:) = glasso_v2(   S_train,L(i),Theta0,fixthese);
        case 'glasso_org'
            Theta(i,:,:) = glasso(      S_train,L(i),Theta0,fixthese);
    end
    % model cardinality
    n0(i) = getZeros(squeeze(   Theta(i,:,:)));
    %n0(i) = sum(sum(sum(abs(Theta(i,:,:))<eps)));
end

Theta = squeeze(Theta);