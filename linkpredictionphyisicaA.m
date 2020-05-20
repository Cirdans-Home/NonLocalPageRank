%% LinkPrediction Problem
% We try to reproduce the results on link-prediction with PageRank from:
% Lü, Linyuan, and Tao Zhou. "Link prediction in complex networks:
% A survey." Physica A: statistical mechanics and its applications 390.6
% (2011): 1150-1170.
clear; clc; close all;

fid = fopen('PhysicaAexample.txt','w+');

problem = ["TestGraphs/gre_216a.mat",...
    "TestGraphs/USAir97.mat",...
    "TestGraphs/netscience.mat",...
    "TestGraphs/power.mat",...
    "TestGraphs/yeast.mat",...
    "TestGraphs/celegansneural"];
for testcase = problem
    % We load the Graph
    load(testcase);
    A = spones(Problem.A);
    flag_of_A_simmetry = issymmetric(A);
    if(flag_of_A_simmetry)
        G = graph(A,'omitselfloops');
    else
        G = digraph(A,'omitselfloops');
    end
    number_of_edges = G.numedges();
    number_of_nodes = G.numnodes();
    e = ones(number_of_nodes,1);
    % We use the K-fold cross-validation with
    rng('default');
    K = 50;
    tau = 0.1;      % We delete 10% of the edges, i.e., the training set
                    % contains 90% of the known links.
    linktoguess = 10;
    
    fprintf('Name of the Network: %s\n',Problem.name);
    fprintf(fid,'Name of the Network: %s\n',Problem.name);
    fprintf(fid,'The network is symmetric : %d ( > 0 is true)\n',flag_of_A_simmetry);
    %% PageRank
    score = zeros(K,1);
    for c = linspace(0.1,0.9,10)
        for k=1:K
            % We pick the 10% of edges at random between all of them
            ind_deleted_edges = randi([1,number_of_edges],...
                floor(tau*number_of_edges),1);
            % We remove the edges we have picked from the graph G
            H = G.rmedge(ind_deleted_edges);
            A = H.adjacency();   % Adjacency matrix of G with removed edges
            % We compute the similarity matrix
            D = 1./(A*e);
            D(D == inf) = 0;
            D = spdiags(D,0,number_of_nodes,number_of_nodes);
            P = D*A;                   % PageRank Initial transition matrix
            % Remove possible zero rows:
            P(P*ones(number_of_nodes,1) == 0,:) = ...
                ones(sum(P*ones(number_of_nodes,1)==0),number_of_nodes)...
                ./number_of_nodes;
            % Compute the similarity matrix by doing (1-c) (I - c*P)^{-1}
            I = eye(number_of_nodes,number_of_nodes);
            X = (1-c)*((I-c*P')\I);
            X = X + X';
            % We do not want to guess the edges we alreay know abot, so we 
            % put their predicted value to -infinity
            X(I>0) = -Inf; X(A>0) = -Inf;
            % Now X contain a ranking of all the edges of the graph, since 
            % X is symmetric we can look at just one of the two triangle 
            % for the rankings
            Xmax = triu(X);
            [max_edges,indmax] = maxk(Xmax(:),linktoguess);
            % We need to go back from the linearized index to the 2D index:
            [I,J] = ind2sub([number_of_nodes number_of_nodes],indmax);
            % The index we want to add are:
            added = [I,J];
            % We check the partialscore 
            partialscore = ...
                ismember(added,G.Edges.EndNodes(ind_deleted_edges,:),'rows');
            score(k) = sum(partialscore)/linktoguess*100;
        end
        
        fprintf(fid,'Standard PageRank c = %1.2f Mean Score = %1.2f +/- %1.2f Median = %1.2f\n',c,mean(score),std(score),median(score));
    end
end
% %% NonLocalPageRank
% alpha = 0.4;
% score = zeros(K,1);
% for k=1:K
%     ind_deleted_edges = folds{k};
%     c = 0.6;
%     % We remove the edges from the graph G
%     H = G.rmedge(ind_deleted_edges);
%     A = G.adjacency();
%     W = H.distances();
%     n = size(A,1);
%     % We compute the similarity matrix
%     W = 1./(W.^alpha);
%     W(W == inf) = 0;
%     D = (1./(W*e));
%     D(D == inf) = 0;
%     D = spdiags(D,0,n,n);
%     P = D*W;
%     % Remove possible zero rows:
%     P(diag(D) == 0,:) = ones(sum(diag(D)==0),n)./n;
%     % Compute the similarity matrix by doing (1-c) (I - c*P)^{-1}
%     I = eye(n,n);
%     X = (1-c)*( (I-c*P.')\I );
%     X = X + X';
%     % We do not want to guess the edges we alreay know of
%     X(I>0) = -Inf; X(A>0) = -Inf;
%     % Now X contain a ranking of all the edges of the graph
%     Xmax = triu(X);
%     [max_edges,indmax] = maxk(Xmax(:),length(ind_deleted_edges));
%     [I,J] = ind2sub([n n],indmax);
%     % The index we want to add are:
%     added = [I,J];
%     partialscore = ismember(added,G.Edges.EndNodes(ind_deleted_edges,:),'rows');
%     score(k) = sum(partialscore)/length(ind_deleted_edges)*100;
% end
%
% fprintf('NonLocal PageRank Mean Score = %1.2f +/- %1.2f\n',mean(score),std(score));