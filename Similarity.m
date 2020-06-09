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
end
if symmetry==true
    X=X+X';
end
end

