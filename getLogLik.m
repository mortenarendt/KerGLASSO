function loglik = getLogLik(Theta,S_test)
%%
%check if Theta is one or several models
if length(size(Theta))>2
    m = size(Theta,1);
else
    m = 1;
end
for i=1:m
    th = squeeze(Theta(i,:,:));
    loglik(i) = log(det(th)) - trace(S_test*th);% + l*log(2*pi);
    if det(th)<0
        'saa er der balada'
    end
end

 