function [tau] = ergodicity(A)
%ERGODICITY Brute-force computation of the ergodicity coefficient of a row
%stochastic matrix
%   A row stochastic matrix
%   tau norm-1 ergodic coefficient of A

n = size(A,1);
tauvec = zeros(n,1);
I = speye(n,n);
e = ones(1,n);

for i=1:n
    tauvec(i) = 0.5*norm(A'*(I - I(:,i)*e),1);
end
tau = max(tauvec);
end

