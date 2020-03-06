%%%%%%%%  Choose penalty based on number of edges

function Lresult = Penalty_based_on_numberofedges_v2(S, Theta0,NfinalEdges,Lstart,fixthese,max_iterations)
p=size(S,1);
UpperBound =(p*p-size(diag(S),1))/2;
if size(Lstart,2)>1
    m= (Lstart(2)-Lstart(1))/2 + Lstart(1);
    Lstart2 = [Lstart(1) m Lstart(2)];
else
    Lstart2=Lstart;
end
%%%% loop starts
Lupdate=Lstart2;
for round=1:max_iterations
    round;
    clear ThetaU n0U TF newmax newmin
    [ThetaU n0U] = glasso_path(S,Lupdate,Theta0,fixthese,'glasso_v3');
    
    % find new penalty space
    n1 = UpperBound - n0U/2;
    if sum(ismember(n1,NfinalEdges))>0 % look at results - break
        break;
    end
    % find new penalty space
    TF = zeros(3,2);
    for i=1:size(n1,2)
        if n1(i)< NfinalEdges;
            TF(i,1)=1;
        elseif  n1(i)> NfinalEdges
            TF(i,2)=1;
        end
    end
    TF = logical(TF);
    newmax = min(Lupdate(TF(:,1)));
    newmin = max(max(Lupdate(TF(:,2))),0.00000001);
    if isempty(newmax)==1;
        newmax = Lupdate(end);
    end
    if isempty(newmin)==1
        newmin = Lupdate(1);
    end
    if newmin==Lupdate(1) & newmax==Lupdate(end);
        break;
    end
    if abs(newmax-newmin) < 0.0001
        break;
    end
    
    
    %Lold = Lupdate;
    %Lupdate = [logspace(log10(newmin),log10(newmax),9)];
    mm = (newmax-newmin)/2 + newmin;
    Lupdate = [newmin mm newmax];
    %nOld = n1;
    
end
idResult = n1==NfinalEdges;
%optimal penalty given NfinalEdges
Lresult = Lupdate(idResult);

if isempty(Lresult)==1;
    Lresult = [newmin newmax];
end
Lresult = min(Lresult(:));
%ThetaResult = squeeze(ThetaU(idResult,:,:))

%G = graph(ThetaResult,'omitselfloops');
%plot(G)


%test
%ThetaR = glasso_path(S,Lresult,Theta0,fixthese,'glasso_v2');

end







