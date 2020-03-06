function [Theta,W,round] = glasso_v3(S,rho,Theta0, fixthese)


fixtheselog = false(size(S,1),1);
fixtheselog(fixthese) = true;
notfixthese = find(~fixtheselog);

[~,p] = size(S);
max_iterations = 1000;% changes from 100 to 1000
t = 1e-4;
convergence_value = t * meanabs(S - diag(diag(S)));

% initialise
%W = S;
W = diag(diag(S + rho*eye(p)));
Theta_old = inv(W);
Theta_old(fixthese,fixthese) = Theta0(fixthese,fixthese);
Theta = Theta_old;
W = inv(Theta_old); 

%%
for round=1:max_iterations
    round
    for j = notfixthese'
        
        i = j;
        
        noti = true(p,1);
        noti(i) = false;
        
        fixtheselog11 = fixtheselog;
        fixtheselog11(i) = [];
        
        % %% Using the DP-GLASSO (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4281944/pdf/nihms637369.pdf)
        % inputs H = Theta11, f = S12*Theta11
        
        % 1) Take out inputs from covariance matrix (W) and precision
        % matrix (Theta)
        w22 = W(i,j);
        theta11 = Theta;
        theta11(i,:) = [];
        theta11(:,j) = [];
        
        s12 = S(:,j);
        s12(i,:) = [];
        
        % 2) Solve Dual via box constraint.
        %X = quadprog(H,f,A,b)
        
        g = quadprog(theta11, s12'*theta11,[],[],[],[], -ones(length(s12),1)*rho,ones(length(s12),1)*rho);
        theta12 = - theta11*(s12 + g) / w22;
        
        theta22 = (1 - (s12 + g)'*theta12) / w22;
        
        w12 = s12 + g;
        
        % 3) reassemble Theta and W
        % 3.1 Theta
        Theta(j,j) = theta22;
        %Theta(j,j) = 1;
        Theta(j,noti) = theta12;
        Theta(noti,j) = theta12;
        % 3.2 W
        W(j,noti) = w12;
        W(noti,j) = w12;
    end
    
    %meanabs(Theta - Theta_old)
    if ~isreal(Theta)
        'problems'
        Theta = NaN(size(Theta));
    end
    
    % Set bad solutions to NA
    %if sum(diag(Theta)) < 0.0001
    %    Theta = NaN(size(Theta));
    %end
    
    %        crit = (meanabs(Theta - Theta_old))
    if (meanabs(Theta - Theta_old) < convergence_value) || sum(sum(isnan(Theta))) >1 %(norm(Theta,1)<t) %Mean of absolute elements of matrix % new stop criteria
        break;
    end
    Theta_old = Theta;
    
end
end

%%
function b = chenLasso(X, Y, lambda, maxIt, tol)
% an algorithm to solve the lasso problem
% from http://pages.stat.wisc.edu/~mchung/teaching/768/matlab/CS/graphicalLasso.m
if nargin < 4, tol = 1e-6; end %nargin = Number of function input arguments
if nargin < 3, maxIt = 1e4; end

% Initialization
[n,p] = size(X);
if p > n
    b = zeros(p,1); % From the null model, if p > n
else
    b = X \ Y;  % From the OLS estimate, if p <= n
end
b_old = b;
i = 0;

% Precompute X'X and X'Y
XTX = X'*X;
XTY = X'*Y;
DD = [];
%XTX = XTX';
% Shooting loop
%set = 1:p;
tr = true(1,p);
while i < maxIt
    i = i+1;
    for j = 1:p
        %jminus = setdiff(1:p,j);
        %jminus = set(tr);
        %S0 = XTX(j,jminus)*b(jminus) - XTY(j);  % S0 = X(:,j)'*(X(:,jminus)*b(jminus)-Y)
        
        tr(j)=0;
        S0 = XTX(j,tr)*b(tr) - XTY(j);
        if S0 > lambda
            b(j) = (lambda-S0) / norm(X(:,j),2)^2;
        elseif S0 < -lambda
            b(j) = -(lambda+S0) / norm(X(:,j),2)^2;
        else
            b(j) = 0;
        end
        tr(j)=1;
    end
    
    delta = norm(b-b_old,1);    % Norm change during successive iterations
    DD = [DD delta];
    %if delta < tol | norm(b,1)<tol, break; end %virker ikke
    if delta < tol || norm(b,1)<0.1, break; end %virker ikke
    b_old = b;
end

if i == maxIt
    fprintf('%s\n', 'Maximum number of iteration reached, shooting may not converge.');
end
end

%%
function A = mysqrtm(B)
%    [V, D] = eig(B);
%    d = diag(D);
%    d(d<0) = 0;
%    A = V*diag(sqrt(d))*V';
[Ut,St,Vt] = svd(B);
A = Ut*sqrt(St)*Vt';
end