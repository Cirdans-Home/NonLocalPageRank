%% LOCALIZATION PHENOMENA FOR THE PAGERANK ALGORITHM

clear; clc;

%% Parameters defininig the cycle
n = 300;         % Number of nodes in the cycle
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
clf
%% Nonlocal pagerank for the graph with added edge
theta = linspace(0,2*pi,n).';
for alpha = [0.01,0.1,0.5,1,100]
    W = Gadd.distances();
    W = 1./(W.^alpha);
    W(W == inf) = 0;
    D = (1./(W*e));
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    P = D*W;
    xGadd = (I - c*P')\((1-c)*e);
    figure(1);
    hold on
    plot(2 +1*cos(theta),2 + xGadd.*sin(theta),'DisplayName',sprintf('\\alpha = %1.2f',alpha),'LineWidth',2);
    hold off
end
hold on
plot(2 +1*cos(theta),2 + 1.*sin(theta),'k--','DisplayName',sprintf('\\alpha = \\infty'),'LineWidth',2);
plot(2 +1*cos(theta),2 + n*xGadd_pr.*sin(theta),'k-.','DisplayName',sprintf('\\alpha = \\infty'),'LineWidth',2);
hold off
set(gca,'FontSize',16)
legend('location','EastOutside');
xticks([]);
yticks([]);
axis tight
axis equal
axes('position',[.25 .30 .35 .40])
box on
for alpha = [0.01,0.1,0.5,1,100]
    W = Gadd.distances();
    W = 1./(W.^alpha);
    W(W == inf) = 0;
    D = (1./(W*e));
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    P = D*W;
    xGadd = (I - c*P')\((1-c)*e);
    figure(1);
    hold on
    plot(2 +1*cos(theta([20:L+7])),2 + xGadd([20:L+7]).*sin(theta([20:L+7])),'LineWidth',2);
    axis tight
    axis equal
    hold off
end
hold on
plot(2 +1*cos(theta([20:L+7])),2 + 1.*sin(theta([20:L+7])),'k--','LineWidth',2);
plot(2 +1*cos(theta([20:L+7])),2 + n*xGadd_pr([20:L+7]).*sin(theta([20:L+7])),'k-.','LineWidth',2);
hold off
xticks([]);
yticks([]);


set(gcf,'color','white')
% export_fig 'cyclegraph.pdf'