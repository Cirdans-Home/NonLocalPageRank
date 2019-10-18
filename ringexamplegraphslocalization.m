%% LOCALIZATION PHENOMENA FOR THE PAGERANK ALGORITHM

clear; clc; close all;

%% Parameters defininig the cycle
n = 100;         % Number of nodes in the cycle
L = 40;          % Node with the added link
c = 0.85;       % Teleportation parameter

% Circulant adjacency matrix:
c0 = [0,1,zeros(1,n-3),1];
A = sparse(gallery('circul',c0));
I = speye(n,n);
e = ones(n,1);

G = digraph(A);                     % Base graph
Gadd = addedge(G,L,1,1);            % Graph with added direct edge
Gadd.Edges.Weight = [];             % Remove the weights

% Standard pagerank
D = 1./(Gadd.adjacency()*e);
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*Gadd.adjacency();
xGadd_pr = powerpr(P,c);

figure(1)
subplot(1,2,1);
%% Nonlocal pagerank for the graph with added edge
theta = linspace(0,2*pi,n).';
for alpha = [0.01,0.1,0.5,1,10]
    W = Gadd.distances();
    W = 1./(W.^alpha);
    W(W == inf) = 0;
    D = (1./(W*e));
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    P = D*W;
    xGadd = (I - c*P.')\((1-c)*e);
    xGadd = xGadd./norm(xGadd,1);
    figure(1);
    subplot(1,2,1);
    hold on
    plot(1:n,xGadd,'DisplayName',sprintf('\\alpha = %1.2f',alpha),'LineWidth',2);
    hold off
end
hold on
plot(1:n,xGadd_pr,'k-.','DisplayName',sprintf('\\alpha = \\infty'),'LineWidth',2);
hold off
legend()
xlabel('Nodes');
ylabel('PageRank');
set(gca,'FontSize',15)

subplot(1,2,2);
h = plot(Gadd,'Layout','circle','LineWidth',1,'MarkerSize',4,'NodeLabel',{});
labelnode(h,[1 10:10:n],[1 10:10:n])
highlight(h,[1 40],'NodeColor','red');
xcoord = (h.XData).';
ycoord = (h.YData).';
axis square tight


