function [X] = Similarity(P,method,c,symmetry,varargin)
%Input: P is a row stochastic matrix
%       Method is :'inverse',
%       c is the PageRank teleportation coefficient
%       If simmetry==True returns symmetric similarity 
%Return Similarity Matrix
[n,~]=size(P);
if strcmp(method,'inverse')
    I=speye(n);
    X=(I-c*P')\I;
    X=(1-c)*X;
elseif strcmp(method,'power')
    X = speye(n);
    for i=1:n
        X(:,i) = powerpr(P,c,X(:,i),1e-3,100,false);
    end
else
    error('Similarity method unkwnown');
end
if symmetry==true
    X=X+X';
end
end

