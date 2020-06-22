% PlotKendal.m - Produces Figure~2 in the Paper from the accompanying
% dataset contained in the matfile kendal.mat
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear all; clc; close all;
load('kendall.mat');
prsize     = size(corr_PRExponential,2);
alphavsize = size(corr_PRExponential,1);
for pp=1:prsize
    figure(1)
    subplot(1,prsize,pp);
    loglog(alphavalue,corr_PRPowerLaw(:,pp)','DisplayName', '$f_{\alpha}(x)= \frac{1}{x^\alpha}$','MarkerSize',7,'LineWidth',2.8);
    hold on
    loglog(alphavalue,corr_PRExponential(:,pp)','DisplayName', '$f_{\alpha}(x)= e^{-\alpha x}$','MarkerSize',7,'LineWidth',2.8);
    title(matrices{pp},'Interpreter','tex','FontSize', 22);
    xlabel('\alpha', 'FontSize', 32);
    ylabel('Kendall \tau', 'FontSize', 32);
    set(gca, ...
        'Box'         , 'on'          , ...
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