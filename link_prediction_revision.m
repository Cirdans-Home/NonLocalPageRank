%% linkprediction.m Tester for the Nonlocal PageRank linkprediction
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear; clc;

%% Loading the dataset
%% Input: Adjacency Matrix Of the graph
datasetchoiche = input(['Select the dataset:\n1) adjnoun\n2)USAir97\ndataset = ']);

switch datasetchoiche
    case 1
        load('TestGraphs/adjnoun.mat');
        dataset = 'adjnoun';
        A=spones(Problem.A);
    case 2
        load('TestGraphs/USAir97.mat')
        %load('TestGraphs/Stranke94.mat')
        dataset = 'USAir97';
        A=spones(Problem.A);
    otherwise
        error('Chose a dataset between 1 and 2');
end

%% Defining the graph and Check If we are considering directed or not directed graph
SimmetryCheck=issymmetric(A);

if SimmetryCheck == 1
   G=graph(A);
else
   G=digraph(A);
end

G = max_connected_subgraph(G);

% %% Plotting the initial Graph
% figure(1)
% subplot(1,3,1);
% initialgraph = plot(G,'layout','force','UseGravity','on');
% xticks([]);
% yticks([]);
% title('Original Graph');

%% Compute Graph Related Quantities
m = numedges(G);
n = numnodes(G);
symmetry=true;
method='inverse';
NumberOfSamples = 30;
alpha_array = [.1 .2 .3 .4 .5 1 2 3 4 5 Inf];  % decay nonlocality
%c_array = [0.08 0.2 0.3 0.4 0.85 0.9]; % pagerank teleportation coeff
c_array = linspace(0.5,0.9,6);
mean_precision=zeros(length(c_array),length(alpha_array));

tau   =  0.1;                       % percentage of removed edges
sigma =  0.5;                       % percentage of edges to predict

for k = 1:length(c_array)
    score = zeros(length(alpha_array),NumberOfSamples);
    for j = 1:NumberOfSamples
        c = c_array(k);                 % pagerank teleportation coeff
        
        ind_deleted_edges = randi([1,m],floor(tau*m),1);                          
        K = sigma*length(ind_deleted_edges);
        
        for i = 1 : length(alpha_array)
            alpha = alpha_array(i);     % decay nonlocality
            [~,score(i,j)] = Predict(G,ind_deleted_edges,c,alpha,K,method,symmetry);
            score(i,j)=score(i,j)/K;
        end
        
        
        if (mod(j,10)==0 || j==1), fprintf('Trial number %d\n', j); end
    end
    
    mean_precision(k,:)=sum(score,2)/NumberOfSamples;
end


h =heatmap(alpha_array,c_array,mean_precision);
h.Title = ['Mean Precision, NoS: ', num2str(NumberOfSamples),' ', dataset,...
           ' tau:',num2str(tau), ' sigma', num2str(sigma)];
h.XLabel = 'NonLocality';
h.YLabel = 'Teleportation';

