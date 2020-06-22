clear all; clc; close all;
load('kendall.mat');
prsize     =size(corr_PRExponential,2);
alphavsize =size(corr_PRExponential,1);
for pp=1:prsize
figure(1)
subplot(1,prsize,pp);
loglog(alphavalue,corr_PRPowerLaw(:,pp)','DisplayName', '$f_{\alpha}(x)= \frac{1}{x^\alpha}$','MarkerSize',7,'LineWidth',2.8);
hold on
loglog(alphavalue,corr_PRExponential(:,pp)','DisplayName', '$f_{\alpha}(x)= e^{-\alpha x}$','MarkerSize',7,'LineWidth',2.8);
%matrices{pp}=erase(matrices{pp},'.mat');
title(matrices{pp},'Interpreter','tex','FontSize', 22);
xlabel('\alpha', 'FontSize', 32);
ylabel('Kendall \tau', 'FontSize', 32); 
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
Lgnd = legend('show');   
set(Lgnd,'Interpreter','latex','Location','northwest')  
Lgnd.FontSize = 20;
hold on
end
hold off
%fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
%Lgnd = legend('show');
%Lgnd.FontSize = 34;
%Lgnd.Position(1) = 0.7;
%Lgnd.Position(2) = 0.001;
set(gcf,'color','white')
tightfig
pause
filename = 'kendall';
savefig(filename);
set(gcf, 'PaperPositionMode', 'auto');
print(filename,'-depsc2');