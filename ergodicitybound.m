%% ERGODICITY WITH RESPECT TO ALPHA
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear; clc;

indexproblem = 1;
alphavalue = logspace(-2,2,40);
typeoftransition = 'powerlaw';   % Input here the type of transformation: exponential or powerlaw
tests = {'USAir97','EuairComplete','barcelona','adjnoun','gre_115','cage9','delaunay_n10','delaunay_n12','3elt'};

for testname = tests
    load(sprintf('TestGraphs/%s.mat',testname{1}));
    A = Problem.A;
    A = A > 0;
    N = size(A,1);
    if issymmetric(A)
        G = graph(A,'omitselfloops');
    else
        G = digraph(A,'omitselfloops');
    end
    
    %% Auxiliary quantities
    e = ones(N,1);
    %% STANDARD PAGERANK
    D = 1./(A*e);
    D(D == inf) = 0;
    Pinf = spdiags(D,0,N,N)*A;
    Ginf = 0.85*Pinf + (1-0.85)/N*(e*e.');
    erginf = ergodicity(Ginf);
    
    %% NONLOCAL PAGERANK
    nodes = (1:N).';
    indexvalue = 1;
    distancematrix = G.distances();
    for alpha = alphavalue
        switch typeoftransition
            case 'powerlaw'
                W = distancematrix;
                W = 1./(W.^alpha);
                W(W == inf) = 0;
                D = (1./(W*e));
                D(D == inf) = 0;
                D = spdiags(D,0,N,N);
                P = D*W;
            case 'exponential'
                W = distancematrix;
                W = exp(-alpha*W);
                W = sparse(W);
                W = W - spdiags(spdiags(W,0),0,N,N);
                D = (1./(W*e));
                D(D == inf) = 0;
                D = spdiags(D,0,N,N);
                P = D*W;
        end
        Galpha = 0.85*P + (1-0.85)/N*(e*e.');
        ergodicityvalue(indexproblem,indexvalue) = erginf - ergodicity(Galpha);
        clear W P Galpha
        indexvalue = indexvalue+1;
    end
    indexproblem = indexproblem+1;
end

%% PLOT RESULTS
figure(1)
heatmap(ergodicityvalue,"XData",alphavalue,"YData",tests);