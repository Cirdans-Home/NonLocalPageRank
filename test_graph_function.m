clear all;
clc;
load('TestGraphs/Stranke94.mat')
A=spones(Problem.A);
SymCheck=issymmetric(A);
%It is an undirected graph
G=graph(A);
m = numedges(G);
n = numnodes(G);
tau=0.1;
ind_deleted_edges = randi([1,m],floor(tau*m),1);
H = G.rmedge(ind_deleted_edges);
AR=adjacency(H);
issymmetric(AR);
%It removes symmetrically the edges!