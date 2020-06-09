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


A = G.adjacency();
G = digraph(spones(A));
G.Edges.Weight = [];
G = max_connected_subgraph(G); % We restrict to the maximum connected
% subgraph

%% Plotting the initial Graph
figure(1)
subplot(1,3,1);
initialgraph = plot(G,'layout','force','UseGravity','on');
xticks([]);
yticks([]);

title('Original Graph');

%% Compute the score for different choices
m = numedges(G);
n = numnodes(G);
NumberOfSamples = 30;
alpha_array = [.1 .2 .3 .4 .5 1 2 3 4 5 Inf];  % decay nonlocality
c_array = [0.01 0.05 0.1 0.4 0.85]; % pagerank teleportation coeff

for k = 1:length(c_array)
    
  
    score = zeros(length(alpha_array),NumberOfSamples);
    for j = 1:NumberOfSamples
        tau = 0.1;                      % percentage of removed edges
        sigma = 1;                      % percentage of edges to predict
        c = c_array(k);                 % pagerank teleportation coeff
        
        ind_deleted_edges = randi([1,m],floor(tau*m),1);                          
        
        for i = 1 : length(alpha_array)
            alpha = alpha_array(i);     % decay nonlocality
            [score(i,j),added] =...
                testingheuristic(G,ind_deleted_edges,alpha,c,sigma);
        end
        
        
        if (mod(j,10)==0 || j==1), fprintf('Trial number %d\n', j); end
    end
    
    numberofdeleted = length(ind_deleted_edges);
    numberofadded = floor(sigma*length(ind_deleted_edges));
    results = score(1:end-1, :)./score(end,:);
    
    % Box-Plot with the results
    figure(2)
    subplot(5,1,k)
    boxplot(results.','Labels',alpha_array(1:end-1));
    ylabel(sprintf('c = %1.2f',c));
    yline(1,'k--');
    if k == 1
        title('Improvement factor with respect to standard rooted PageRank')
    elseif k == length(c_array)
        xlabel('alpha')
    end
end

%% Plot the resulting graph
H = G.rmedge(ind_deleted_edges);
Gnew = H.addedge(added(:,1),added(:,2));
figure(1)
subplot(1,3,2);
plot(G,'XData',initialgraph.XData,'YData',initialgraph.YData)
xticks([]);
yticks([]);
title('Graph with removed edges');
subplot(1,3,3);
newgraph = plot(Gnew,'XData',initialgraph.XData,'YData',...
    initialgraph.YData);
highlight(newgraph,added(:,1),added(:,2),'EdgeColor','r',...
    'LineWidth',2);
xticks([]);
yticks([]);
title('Graph with added edges')

