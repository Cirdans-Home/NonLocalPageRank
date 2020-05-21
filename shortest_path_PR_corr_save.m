%% TESTING DISTANCES AND SMOOTHING FUNCTIONs
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute
% Revision 21th May 2020
clear; clc; close all;
matrices={'USAir97.mat','netscience.mat','bitcoin.mat'};
%'pesa.mat'
alphavalue=logspace(-3,1,20);
corr_PRPowerLaw   =zeros(size(alphavalue,2),length(matrices)); 
corr_PRExponential=zeros(size(alphavalue,2),length(matrices));

for pp=1:length(matrices)
load(['TestGraphs/',matrices{pp}]);
A = spones(Problem.A);
N = size(A,1);
A=A-spdiags(diag(A),0,N,N);
flag_of_A_simmetry = issymmetric(A);
if(flag_of_A_simmetry)
        G = graph(A,'omitselfloops');
else
        G = digraph(A,'omitselfloops');
end
e = ones(N,1);
%% Shortest Path Distance
W = G.distances();
%% STANDARD PAGERANK
D = 1./(A*e);
D(D == inf) = 0;
Pinf = spdiags(D,0,N,N)*A;
[PageRank,flag,reshist] = powerpr(Pinf);
%% NONLOCAL PAGERANK
% Exponential
PagerankPowerLawSPD=zeros(N,size(alphavalue,2));
for i=1:size(alphavalue,2)
            alpha = alphavalue(i);
            W1 = 1./(W.^alpha);
            W1(W1 == inf) = 0;
            D1 = (1./(W1*e));
            D1(D1 == inf) = 0;
            D1 = spdiags(D1,0,N,N);
            P1 = D1*W1;
            [PagerankPowerLawSPD(:,i),flag,reshist] = powerpr(P1);
end
PagerankExponenetialSPD=zeros(N,size(alphavalue,2));
for i=1:size(alphavalue,2)
            alpha = alphavalue(i);
            W1 = exp(-alpha*W);
            W1 = sparse(W1);
            W1 = W1 - spdiags(spdiags(W1,0),0,N,N);
            D1 = (1./(W1*e));
            D1(D1 == inf) = 0;
            D1 = spdiags(D1,0,N,N);
            P1 = D1*W1;
            [PagerankExponenetialSPD(:,i),flag,reshist] = powerpr(P1);
end

for i=1:size(alphavalue,2)
    corr_PRPowerLaw(i,pp)   =corr(PageRank,PagerankPowerLawSPD(:,i), 'Type','Kendall');
    corr_PRExponential(i,pp)=corr(PageRank,PagerankExponenetialSPD(:,i), 'Type','Kendall');
end
end
filename = 'kendall';
save(filename,'corr_PRPowerLaw', 'corr_PRExponential','alphavalue','matrices');