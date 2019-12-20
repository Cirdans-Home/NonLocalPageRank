%% TESTING DISTANCES AND SMOOTHING FUNCTIOND
clear; clc;
load('TestGraphs/USAir97.mat');
A = Problem.A;
A = A > 0;
A=spones(A);
N = size(A,1);
G = digraph(A,'omitselfloops');
%% Logarithmic distance Matrix
e = ones(N,1);
I = speye(N,N);
L = spdiags(A*e,0,N,N) - A;
S = inv(I+L);
[row,col] = find(S < 0); %Check if there are numerical negative entries
S(row,col) = 0;
H = log(S);
h = diag(H);
U = h*e' - H;
WLog = 0.5*(U+U.');
%% Shortest Path Distance
W = G.distances();
%% STANDARD PAGERANK
D = 1./(A*e);
D(D == inf) = 0;
Pinf = spdiags(D,0,N,N)*A;
[PageRank,flag,reshist] = powerpr(Pinf);
%% NONLOCAL PAGERANK
alphavalue=[0.1,0.5,1,3];
% Exponential
PagerankPowerLawSPD=zeros(N,size(alphavalue,2));
PagerankPowerLawLOG=zeros(N,size(alphavalue,2));
for i=1:size(alphavalue,2)
            alpha = alphavalue(i);
            W1 = 1./(W.^alpha);
            W1(W1 == inf) = 0;
            D1 = (1./(W1*e));
            D1(D1 == inf) = 0;
            D1 = spdiags(D1,0,N,N);
            P1 = D1*W1;
            [PagerankPowerLawSPD(:,i),flag,reshist] = powerpr(P1);
            
            W2 = 1./(WLog.^alpha);
            W2(W2 == inf) = 0;
            D2 = (1./(W2*e));
            D2(D2 == inf) = 0;
            D2 = spdiags(D2,0,N,N);
            P2 = D2*W2;
            [PagerankPowerLawLOG(:,i),flag,reshist] = powerpr(P2);
            
end
PagerankExponenetialSPD=zeros(N,size(alphavalue,2));
PagerankExponenetialLOG=zeros(N,size(alphavalue,2));
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
            
            W2 = exp(-alpha*WLog);
            W2 = sparse(W2);
            W2 = W2 - spdiags(spdiags(W2,0),0,N,N);
            D2 = (1./(W2*e));
            D2(D2 == inf) = 0;
            D2 = spdiags(D2,0,N,N);
            P2 = D2*W2;
            [PagerankExponenetialLOG(:,i),flag,reshist] = powerpr(P2);
end
PagerankPowerLawSPD1=zeros(N,size(alphavalue,2));
PagerankPowerLawLOG1=zeros(N,size(alphavalue,2));
PagerankExponenetialSPD1=zeros(N,size(alphavalue,2));
PagerankExponenetialLOG1=zeros(N,size(alphavalue,2));
for i=1:size(alphavalue,2)
    PagerankPowerLawSPD1(:,i)=linearize(PagerankPowerLawSPD(:,i),0,1);
    PagerankPowerLawLOG1(:,i)=linearize(PagerankPowerLawLOG(:,i),0,1);
    PagerankExponenetialSPD1(:,i)=linearize(PagerankExponenetialSPD(:,i),0,1);
    PagerankExponenetialLOG1(:,i)=linearize(PagerankExponenetialLOG(:,i),0,1);
end
PageRank1=linearize(PageRank,0,1);
figure(1);
for i=1:size(alphavalue,2)
    subplot(2,size(alphavalue,2),i);
    loglog(PageRank1,PagerankPowerLawSPD1(:,i),'rx','DisplayName', '$f_{\alpha}(x)= \frac{1}{x^\alpha}$','MarkerSize',7,'LineWidth',2.8);
    set(gca, ...
  'Box'         , 'on'          , ...
  'PlotBoxAspectRatio',[2 2 2]  , ...
  'TickDir'     , 'in'          , ...
  'TickLength'  , [.01 .01]     , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'YGrid'       , 'off'          , ...
  'XGrid'       , 'off'         , ...
  'XColor'      , [.3 .3 .3]    , ...
  'YColor'      , [.3 .3 .3]    , ...
  'LineWidth'   , 1            , ...        
  'FontSize'    , 12);
   title(['\alpha = ',num2str(alphavalue(i))],'Interpreter','tex');
   if i==1
       xlabel('PageRank', 'FontSize', 25);
       ylabel('Nonlocal PageRank', 'FontSize', 25);
   end
end
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('show');
Lgnd.FontSize = 34;
set(Lgnd,'Interpreter','latex')
Lgnd.Position(1) = 0.35;
Lgnd.Position(2) = 0.875;
for i=1:size(alphavalue,2)
    subplot(2,size(alphavalue,2),i+size(alphavalue,2));
    loglog(PageRank1,PagerankExponenetialSPD1(:,i),'bo','DisplayName', '$f_{\alpha}(x)= e^{-\alpha x}$','MarkerSize',7,'LineWidth',2.8);
    set(gca, ...
  'Box'         , 'on'          , ...
  'PlotBoxAspectRatio',[2 2 2]  , ...
  'TickDir'     , 'in'          , ...
  'TickLength'  , [.01 .01]     , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'YGrid'       , 'off'          , ...
  'XGrid'       , 'off'         , ...
  'XColor'      , [.3 .3 .3]    , ...
  'YColor'      , [.3 .3 .3]    , ...
  'LineWidth'   , 1             , ...        
  'FontSize'    , 12);
   title(['\alpha = ',num2str(alphavalue(i))],'Interpreter','tex');
   if i==1
       xlabel('PageRank', 'FontSize', 25);
       ylabel('Nonlocal PageRank', 'FontSize', 25);
   end
end
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('show');
Lgnd.FontSize = 34;
set(Lgnd,'Interpreter','latex')
Lgnd.Position(1) = 0.375;
Lgnd.Position(2) = 0.405;
set(gcf,'color','white')
tightfig
pause
filename = 'shortestpath';
savefig(filename);
set(gcf, 'PaperPositionMode', 'auto');
print(filename,'-depsc2');


figure(2)
for i=1:size(alphavalue,2)
    subplot(2,size(alphavalue,2),i);
    loglog(PagerankPowerLawSPD1(:,i),PagerankPowerLawLOG1(:,i),'rx','DisplayName', '$f_{\alpha}(x)= \frac{1}{x^\alpha}$','MarkerSize',5,'LineWidth',2.8);
    set(gca, ...
  'Box'         , 'on'          , ...
  'PlotBoxAspectRatio',[2 2 2]  , ...
  'TickDir'     , 'in'          , ...
  'TickLength'  , [.01 .01]     , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'YGrid'       , 'off'          , ...
  'XGrid'       , 'off'         , ...
  'XColor'      , [.3 .3 .3]    , ...
  'YColor'      , [.3 .3 .3]    , ...
  'LineWidth'   , 1             , ...        
  'FontSize'    , 12);
   title(['\alpha = ',num2str(alphavalue(i))],'Interpreter','tex');
   if i==1
       xlh=xlabel('Shortest Path', 'FontSize', 25);
       ylabel('Logarithmic', 'FontSize', 25);
   end
end
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.FontSize = 34;
set(Lgnd,'Interpreter','latex')
Lgnd.Position(1) = 0.35;
Lgnd.Position(2) = 0.875;
set(gcf,'color','white')
for i=1:size(alphavalue,2)
    subplot(2,size(alphavalue,2),i+size(alphavalue,2));
    loglog(PagerankExponenetialSPD1(:,i),PagerankExponenetialLOG1(:,i),'bo','DisplayName', '$f_{\alpha}(x)= e^{-\alpha x}$','MarkerSize',5,'LineWidth',2.8);
    set(gca, ...
  'Box'         , 'on'          , ...
  'PlotBoxAspectRatio',[2 2 2]  , ...
  'TickDir'     , 'in'          , ...
  'TickLength'  , [.01 .01]     , ...
  'XMinorTick'  , 'off'         , ...
  'YMinorTick'  , 'off'         , ...
  'YGrid'       , 'off'          , ...
  'XGrid'       , 'off'         , ...
  'XColor'      , [.3 .3 .3]    , ...
  'YColor'      , [.3 .3 .3]    , ...
  'LineWidth'   , 1             , ...        
  'FontSize'    , 12);
   title(['\alpha = ',num2str(alphavalue(i))],'Interpreter','tex');
   if i==1
       xlh=xlabel('Shortest Path', 'FontSize', 25);
       ylabel('Logarithmic', 'FontSize', 25);
   end
end
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('show');
Lgnd.FontSize = 34;
set(Lgnd,'Interpreter','latex')
Lgnd.Position(1) = 0.375;
Lgnd.Position(2) = 0.405;
set(gcf,'color','white')

tightfig
pause
filename = 'exp_vs_log';
savefig(filename);
set(gcf, 'PaperPositionMode', 'auto');
print(filename,'-depsc2');