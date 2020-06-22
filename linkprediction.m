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

TestGraph={'adjnoun','USAir97','cage9', 'delaunay_n12','3elt',...
    'EuairComplete','delaunay_n10'};

Result=zeros(length(TestGraph),2);
Results_TAB = [];
% We initizialize the random number generator to have reproducible results
rng(100)

for problem=1:length(TestGraph)
    
    load(['TestGraphs/',TestGraph{problem},'.mat']);
    A=spones(Problem.A);
    % Defining the graph: 
    % Observe that Every NonSymmetric Graph Will be symmetrized
    SimmetryCheck=issymmetric(A);
    
    if SimmetryCheck==0
        A=A+A';
        A=spones(A);
        SimmetryCheck=issymmetric(A);
    end
    
    if SimmetryCheck == 1
        G=graph(A,'omitselfloops');
    else
        G=digraph(A);
    end
    
    % Compute Graph Related Quantities
    m = numedges(G);
    n = numnodes(G);
    symmetry=true;
    alpha_array = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 2 3 5 Inf];  % decay nonlocality
    c_array = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.85 0.9 0.99 ];  % pagerank teleportation coeff

    tau         =  0.1;  % percentage of removed edges
    fprintf('Dataset = %s \t |V| = %d \t |E| = %d \t |E|/|V| = %1.5f\n\n',TestGraph{problem},n,m,m/(n^2));
    NumTrial1 = 15;
    numfolds = 10;
    LocalScores = zeros(length(TestGraph), NumTrial1);
    NLocalScores= zeros(length(TestGraph),NumTrial1);
    
    
    %%% Variables for final results table
    num_nodes = n;
    num_edges = m;
    density = m/(n^2);
    dataset = string(TestGraph{problem});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for nmt1=1:NumTrial1
        ind_deleted_edges_ToTrain=randperm(m);
        ind_deleted_edges_ToTrain=ind_deleted_edges_ToTrain(1:floor(tau*m));
        KT=length(ind_deleted_edges_ToTrain);
        GT=G.rmedge(ind_deleted_edges_ToTrain);
        mT = numedges(GT);
        nT = numnodes(GT);
        
        arch=[1:mT];
        cV = cvpartition(arch,'KFold',numfolds,'Stratify',false);
        best_c=zeros(cV.NumTestSets,1);
        best_alpha=zeros(cV.NumTestSets,1);
        best_score=zeros(cV.NumTestSets,1);
        best_c_L=zeros(cV.NumTestSets,1);
        best_score_L=zeros(cV.NumTestSets,1);
        for nmt2=1:cV.NumTestSets
            
            ind_deleted_edges = find(cV.test(nmt2)==1);
            K = length(ind_deleted_edges);
            
            score = zeros(length(c_array),length(alpha_array));
            
            H = G.rmedge(ind_deleted_edges);
            AR = H.adjacency(); %Adjacency Matrix after removing edges
            W = H.distances();
            
            for j = 1:length(alpha_array)
                alpha = alpha_array(j);     % decay nonlocality
                
                P = make_matrix_P(W,AR,alpha);
                
                for i = 1:length(c_array)
                    c = c_array(i);
                    [~,score(i,j)] = FT_Predict(GT,AR,P,c,K,symmetry,ind_deleted_edges);
                    score(i,j)=score(i,j)/K;
                end
                
            end
            
            % NonLocal
            M=max(max(score(:,1:end-1)));
            [bi,bj]=find(score(:,1:end-1)==M);
            best_c(nmt2)=c_array(bi(1));
            best_alpha(nmt2)=alpha_array(bj(1));
            best_score(nmt2)=M;
            % Local
            [aL,bL]=max(score(:,end));
            best_c_L(nmt2)=c_array(bL(1));
            best_score_L(nmt2)=aL(1);
            
        end
        
        % NonLocal
        [Snl,INTSnl]=max(best_score);
        cnl=best_c(INTSnl); % best value of c for this trial
        alpha1=best_alpha(INTSnl); % best value of alpha for this trial
        % Local
        [Sl,INTSl]=max(best_score_L);
        cl=best_c_L(INTSl); % best value of c for this trial
        [~,NLocalScores(problem,nmt1)] = Predict(G,ind_deleted_edges_ToTrain,...
            cnl(1),alpha1,KT,'inverse',symmetry);
        NLocalScores(problem,nmt1)=NLocalScores(problem,nmt1)/KT;
        [~,LocalScores(problem,nmt1)] = Predict(G,ind_deleted_edges_ToTrain,...
            cl(1),Inf,KT,'inverse',symmetry);
        LocalScores(problem,nmt1)=LocalScores(problem,nmt1)/KT;
        
        %%% Variables for final results table %%%%%%%%%%
        trial = nmt1;
        accuracy_local = LocalScores(problem,nmt1);
        accuracy_nonlocal = NLocalScores(problem,nmt1);
        alpha = alpha1;
        c_nonlocal = cnl;
        c_local = cl;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Results_TAB = [Results_TAB;
            table(dataset, ...
            num_nodes, ...
            num_edges, ...
            density, ...
            trial, ...
            numfolds, ...
            alpha, ...
            c_nonlocal, ...
            accuracy_nonlocal, ...
            c_local, ...
            accuracy_local)];
        disp(Results_TAB(Results_TAB.dataset == dataset,5:end))
        
        filepath = sprintf('results/results_cv_%s.csv', string(date));
        writetable(Results_TAB,filepath)
    end
    LM=sum(LocalScores(problem,:),2)/NumTrial1;
    NLM=sum(NLocalScores(problem,:),2)/NumTrial1;

    Result(problem,1)=LM;
    Result(problem,2)=NLM;
    
end
save KFoldCVM Results_TAB Result LocalScores NLocalScores
filepath = sprintf('results/results_cv_%s.csv', string(date));
writetable(Results_TAB,filepath)