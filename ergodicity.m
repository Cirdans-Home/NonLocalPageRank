function [tau] = ergodicity(A)
%ERGODICITY Brute-force computation of the ergodicity coefficient of a row
%stochastic matrix
%   A row stochastic matrix
%   tau norm-1 ergodic coefficient of A
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute


n = size(A,1);
% tauvec = zeros(n,1);
I = speye(n,n);
e = ones(1,n);

% for i=1:n
%     tauvec(i) = 0.5*norm(A'*(I - I(:,i)*e),1);
% end
tauvec = arrayfun(@(i) 0.5*norm(A'*(I - I(:,i)*e),1),1:n);
tau = max(tauvec);
end

