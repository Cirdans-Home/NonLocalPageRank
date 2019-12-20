%% ERGODICITY WITH RESPECT TO ALPHA

clear; clc; 

load('TestGraphs/adjnoun.mat');
A = Problem.A;
A = A > 0;
N = size(A,1);
G = graph(A,'omitselfloops');

%% Auxiliary quantities
e = ones(N,1);
%% STANDARD PAGERANK
D = 1./(A*e);
D(D == inf) = 0;
Pinf = spdiags(D,0,N,N)*A;
Ginf = 0.85*Pinf + (1-0.85)/N*(e*e.');
erginf = ergodicity(Ginf);
[PageRank,flag,reshist] = powerpr(Pinf);

%% NONLOCAL PAGERANK
ranks = []; ergodicityvalue = []; ErrPr = []; ErrNor = []; ergodicityvalue_estimate = [];
alphavalue = (0.1:0.1:40);
typeoftransition = 'exponential';   % Input here the type of transformation: exponential or powerlaw
nodes = (1:N).';

for alpha = alphavalue
    switch typeoftransition
        case 'powerlaw'
            W = G.distances();
            W = 1./(W.^alpha);
            W(W == inf) = 0;
            D = (1./(W*e));
            D(D == inf) = 0;
            D = spdiags(D,0,N,N);
            P = D*W;
        case 'exponential'
            W = G.distances();
            W = exp(-alpha*W);
            W = sparse(W);
            W = W - spdiags(spdiags(W,0),0,N,N);
            D = (1./(W*e));
            D(D == inf) = 0;
            D = spdiags(D,0,N,N);
            P = D*W;
    end
    [PageRankalpha,flag,reshist] = powerpr(P);
    PageRankalpha = PageRankalpha./norm(PageRankalpha,1);
    Galpha = 0.85*P + (1-0.85)/N*(e*e.');
    ergodicityvalue = [ergodicityvalue,ergodicity(Galpha)];
    ergodicityvalue_estimate = [ergodicityvalue_estimate,norm(Galpha,1)];
    ErrNor = [ErrNor,norm(Galpha - Ginf,1)];
    ErrPr = [ErrPr,norm(PageRankalpha - PageRank,1)];
    clear W P Galpha
end

%% PLOT RESULTS

figure(1);
subplot(1,2,1);
loglog(alphavalue,1./(1-ergodicityvalue),'b-',...
    alphavalue,1./(1-erginf*ones(size(alphavalue))),'r--','LineWidth',2);
legend({'Non--Local','Local'},'Location','SouthEast')
xlabel('$\alpha$','Interpreter','LaTeX');
ylabel('$\tau_1(G_\alpha)$','Interpreter','LaTeX');
axis([alphavalue(1) alphavalue(end) 1/(1-min(ergodicityvalue)) (1/(1-min(erginf))+1)])
axis square
set(gca,'FontSize',14)

figure(1);
subplot(1,2,2);
loglog(alphavalue,ErrPr,'b-',...
    alphavalue,ErrNor./(1-ergodicityvalue),'r--',...
    'LineWidth',2);
xlabel('$\alpha$','Interpreter','LaTeX');
ylabel('$\|s_\infty - s_\alpha \|_1$','Interpreter','LaTeX');
legend({'Error','Bound'},'Location','SouthWest')
axis tight
axis square
set(gca,'FontSize',14)

set(gcf,'color','white')