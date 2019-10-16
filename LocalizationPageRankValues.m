%% LOCALIZATION OF THE PAGERANK MEASURE
% Examples for the stability and non locality section of the paper

clear; clc;

load('TestGraphs/adjnoun.mat');
A = Problem.A;
G = digraph(A);
n = size(A,1);
e = ones(n,1);

figure(100);
h = plot(G,'Layout','subspace');
xy(:,1) = h.XData;
xy(:,2) = h.YData;

c = 0.90;

% Nonlocal pagerank for alpha = 0.01
alpha = 0.01;
W = G.distances();
W = 1./(W.^alpha);
W(W == inf) = 0;
D = (1./(W*e));
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*W;
x1 = powerpr(P,c);
% Nonlocal pagerank for alpha = 0.5
alpha = 0.5;
W = G.distances();
W = 1./(W.^alpha);
W(W == inf) = 0;
D = (1./(W*e));
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*W;
x2 = powerpr(P,c);
% Nonlocal pagerank for alpha = 3
alpha = 3;
W = G.distances();
W = 1./(W.^alpha);
W(W == inf) = 0;
D = (1./(W*e));
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*W;
x3 = powerpr(P,c);
% Standard pagerank
D = 1./(A*e);
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*A;
x4 = powerpr(P,c);

scalingnodesize = 15*1/max(max([x1,x2,x3,x4]));

figure(1);
subplot(2,2,1)
plot(G,'XData',xy(:,1),'YData',xy(:,2),'NodeCData',x1,'MarkerSize',scalingnodesize*x1,'NodeLabel',{})
axis tight
axis square
title('\alpha = 0.01');
xticks([]);
yticks([]);

subplot(2,2,2)
plot(G,'XData',xy(:,1),'YData',xy(:,2),'NodeCData',x2,'MarkerSize',scalingnodesize*x2,'NodeLabel',{})
axis tight
axis square
title('\alpha = 0.5');
xticks([]);
yticks([]);

subplot(2,2,3)
plot(G,'XData',xy(:,1),'YData',xy(:,2),'NodeCData',x3,'MarkerSize',scalingnodesize*x3,'NodeLabel',{})
axis tight
axis square
title('\alpha = 3');
xticks([]);
yticks([]);


subplot(2,2,4)
plot(G,'XData',xy(:,1),'YData',xy(:,2),'NodeCData',x4,'MarkerSize',scalingnodesize*x4,'NodeLabel',{})
axis tight
axis square
title('\alpha = \infty');
xticks([]);
yticks([]);

set(gcf,'color','white')