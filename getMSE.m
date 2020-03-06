function [MSE MSEp] = getMSE(A,Ahat)

for i=1:size(Ahat,1)
    e = A - squeeze(Ahat(i,:,:));
    e = e(~isnan(e));
    MSE(i) = sum(e.^2) / length(e);
end
e =  A - A(randperm(size(A,1)),:);
e = e(~isnan(e));
MSEp = sum(e.^2) / length(e);
