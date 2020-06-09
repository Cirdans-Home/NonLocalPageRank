%% linkprediction.m Tester for the Nonlocal PageRank linkprediction
% This tester reproduces the performances of the non-local pagerank on the
% problem of predicting the edges of a given graph (digraph) G. Namely, we
% suppose having a snapshot of the graph (digraph) G at time t0, and we
% know that at time t1 a new graph (digraph) G' will need to be computed in
% such way that G and G' have the same nodes, but some edges will have been
% added to G. To test the performances the present code will start from a
% graph (digraph) G, and will remove a chosen amount of edges from it, than
% the heuristic algorithm will try to guess what edges to add. A score
% based on the number of correctedly guessed links is then produced.
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear; clc;

%% Loading the dataset

datasetchoiche = input(['Select the dataset:\n1) Barcelona\n2) USAir97\ndataset = ']);

switch datasetchoiche
    case 1
        load('TestGraphs/barcelona.mat');
        dataset = 'Barcelona';
    case 2
        load('TestGraphs/USAir97.mat')
        dataset = 'USAir97';
        A = spones(Problem.A);
        G = graph(A);
    otherwise
        error('Chose a dataset between 1 and 2');
end
