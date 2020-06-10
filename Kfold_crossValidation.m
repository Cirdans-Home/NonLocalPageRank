%% linkprediction.m Tester for the Nonlocal PageRank linkprediction
%  Kfolding Training UNDIRECTED GRAPH
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear; clc;

%% Loading the dataset

% %% Input: Adjacency Matrix Of the graph
% datasetchoiche = input(['Select the dataset:\n1)adjnoun\n2)USAir97',...
%                          '\n3)EUair\n4)Tube\ndataset =']);
% 
% switch datasetchoiche
%     case 1
%         load('TestGraphs/adjnoun.mat');
%         dataset = 'adjnoun';
%         A=spones(Problem.A);
%     case 2
%         load('TestGraphs/USAir97.mat')
%          dataset = 'USAir97';
%         A=spones(Problem.A);
%     case 3
%         load('TestGraphs/EuairComplete.mat')
%         dataset = 'EuAir';
%         A=spones(Problem.A);
%     case 4
%         load('TestGraphs/TubeComplete.mat')
%         dataset = 'Tube';
%         A=spones(Problem.A);
%     otherwise
%         error('Chose a dataset between 1 and 4');
% end

TestGraph={'adjnoun', 'USAir97', 'EuairComplete','football',...
           'celegans_metabolic','polbooks','mycielskian9',...
           'delaunay_n10','data','USpowerGrid'};
% TestGraph={'adjnoun'};
Result=zeros(length(TestGraph),2);

Results_TAB = [];


rng(100)


for problem=1:length(TestGraph)

load(['TestGraphs/',TestGraph{problem},'.mat']);
A=spones(Problem.A);

%% Defining the graph and Check If we are considering directed or not directed graph
SimmetryCheck=issymmetric(A);

if SimmetryCheck == 1
   G=graph(A,'omitselfloops');
else
   G=digraph(A);
end

%G = max_connected_subgraph(G);

%% Compute Graph Related Quantities
m = numedges(G);
n = numnodes(G);
symmetry=true;
method='inverse';
alpha_array = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 2 3 5 Inf];  % decay nonlocality
c_array = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.85 0.9 0.99 ];  % pagerank teleportation coeff

%alpha_array = [.4 .5 Inf];   
%c_array     = [.1 .85 .9]; 
tau         =  0.1;  % percentage of removed edges
fprintf('Dataset = %s \t |V| = %d \t |E| = %d \t |E|/|V| = %1.5f\n\n',TestGraph{problem},n,m,m/(n^2));
NumTrial1=10;
LocalScores = zeros(length(TestGraph), NumTrial1);
NLocalScores= zeros(length(TestGraph),NumTrial1);

num_nodes = n;
num_edges = m;
density = m/(n^2);
dataset = string(TestGraph{problem});


for nmt1=1:NumTrial1
%     fprintf('Trial n: %d\n\n',nmt1);
    ind_deleted_edges_ToTrain=randperm(m);
    ind_deleted_edges_ToTrain=ind_deleted_edges_ToTrain(1:floor(tau*m));
    KT=length(ind_deleted_edges_ToTrain);
    GT=G.rmedge(ind_deleted_edges_ToTrain);
    mT = numedges(GT);
    nT = numnodes(GT);

    arch=[1:mT];
    cV = cvpartition(arch,'KFold',10,'Stratify',false);
    best_c=zeros(cV.NumTestSets,1);
    best_alpha=zeros(cV.NumTestSets,1);
    best_score=zeros(cV.NumTestSets,1);
    best_c_L=zeros(cV.NumTestSets,1);
    best_score_L=zeros(cV.NumTestSets,1);
    for nmt2=1:cV.NumTestSets

        ind_deleted_edges = find(cV.test(nmt2)==1);                          
        K = length(ind_deleted_edges);
        %Gbar=GT.rmedge(ind_deleted_edges);
        score = zeros(length(c_array),length(alpha_array));
        for i = 1:length(c_array)
        c = c_array(i);                 % pagerank teleportation coeff
            for j = 1 : length(alpha_array)
            alpha = alpha_array(j);     % decay nonlocality
            [~,score(i,j)] = Predict(GT,ind_deleted_edges,c,alpha,K,method,symmetry);
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
                             cnl(1),alpha1,KT,method,symmetry);
    NLocalScores(problem,nmt1)=NLocalScores(nmt1)/KT;
    [~,LocalScores(problem,nmt1)] = Predict(G,ind_deleted_edges_ToTrain,...
                            cl(1),Inf,KT,method,symmetry);
    LocalScores(problem,nmt1)=LocalScores(nmt1)/KT;
    
%     fprintf('alpha  =  %1.3f \t cnl = %1.3f \t cl =%1.3f \n', alpha1, cnl, cl)
%     fprintf('Local = %1.5f \t Nonlocal = %1.5f\n\n', LocalScores(problem,nmt1),  NLocalScores(problem,nmt1));

    trial = nmt1;
    accuracy_local = LocalScores(problem,nmt1);
    accuracy_nonlocal = NLocalScores(problem,nmt1);
    alpha = alpha1;
    c_nonlocal = cnl;
    c_local = cl;
    
    Results_TAB = [Results_TAB; table(dataset, num_nodes, num_edges, density, trial, alpha, c_nonlocal, accuracy_nonlocal, c_local, accuracy_local)];
    disp(Results_TAB(:,5:end))
end
LM=sum(LocalScores(problem,:),2)/NumTrial1;
NLM=sum(NLocalScores(problem,:),2)/NumTrial1;
%h =heatmap(alpha_array,c_array,mean_precision);
%h.Title = ['Mean Precision, NoS: ', num2str(NumberOfSamples),' ', dataset,...
           %' Tau:',num2str(tau), ' Sigma:', num2str(sigma)];
%h.XLabel = 'NonLocality';
%h.YLabel = 'Teleportation';

Result(problem,1)=LM;
Result(problem,2)=NLM;

end
save KFoldCVM Results_TAB Result LocalScores NLocalScores
writetable(Results_TAB,'results.csv')
