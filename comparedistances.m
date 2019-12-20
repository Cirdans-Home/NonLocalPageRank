%% COMPARE DISTANCES FOR A DIGRAPH \Gamma
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute


clear; clc; close all;

load('TestGraphs/gre_115.mat');
A = spones(Problem.A);
n = size(A,1);

% Building the digraph:
G = digraph(A,'omitselfloops');

% Shortest-path distance
W1 = G.distances();

% Logarithmic distance
e = ones(n,1);
I = speye(n,n);
L = spdiags(A*e,0,n,n) - A;
S = inv(I+L);
[row,col] = find(S < 0); %Check if there are numerical negative entries
S(row,col) = 0;
H = log(S);
h = diag(H);
U = h*e' - H;

W2 = 0.5*(U+U.');

%% PLOT ROUTINES

figure(1);
[ha, pos] = tight_subplot(1, 3,[.01 .03],[.1 .01],[.01 .01]);
axes(ha(1)); spy(A,'s'); xticks([]); yticks([]); title('Adjacency'); set(gca,'FontSize',14)
axes(ha(2)); spyc(W1,'autumn',1); xticks([]); yticks([]); title('Shortest Path'); set(gca,'FontSize',14)
axes(ha(3)); spyc(W2,'autumn',1); xticks([]); yticks([]); title('Logarithmic'); set(gca,'FontSize',14)
set(gcf,'color','white');

figure(2);
[hb, posb] = tight_subplot(1, 2,[.01 .03],[.1 .01],[.01 .01]);
axes(hb(1)); spy(L);
axes(hb(2)); spyc(S,'autumn',1); xticks([]); yticks([]); title('Logarithmic'); set(gca,'FontSize',14)
