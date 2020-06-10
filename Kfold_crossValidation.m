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

%TestGraph={'adjnoun', 'USAir97', 'EuairComplete','football',...
           %'celegans_metabolic','polbooks','mycielskian9',...
           %'delaunay_n10','data','USpowerGrid'};
TestGraph={'USAir97'};
Result=zeros(length(TestGraph),2);
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
alpha_array = [.1 .2 .3 .4 .5 1 2 3 4 5 6 7 Inf];  % decay nonlocality
c_array = [0.08 0.2 0.3 0.4 0.5 0.6 0.7 0.85 0.9 0.99]; % pagerank teleportation coeff

%alpha_array = [.4 .5 Inf];   
%c_array     = [.1 .85 .9]; 
tau         =  0.1;  % percentage of removed edges
TestGraph{problem}
NumTrial1=10;
LocalScores = zeros(length(TestGraph), NumTrial1);
NLocalScores= zeros(length(TestGraph),NumTrial1);
for nmt1=1:NumTrial1
    ['Trial n: ',num2str(nmt1)]
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
    for  nmt2=1:cV.NumTestSets

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
    cnl=best_c(INTSnl);
    alpha1=best_alpha(INTSnl);
    % Local
    [Sl,INTSl]=max(best_score_L);
    cl=best_c_L(INTSl);
    [~,NLocalScores(problem,nmt1)] = Predict(G,ind_deleted_edges_ToTrain,...
                             cnl(1),alpha1,KT,method,symmetry);
    NLocalScores(problem,nmt1)=NLocalScores(nmt1)/KT;
    [~,LocalScores(problem,nmt1)] = Predict(G,ind_deleted_edges_ToTrain,...
                            cl(1),Inf,KT,method,symmetry);
    LocalScores(problem,nmt1)=LocalScores(nmt1)/KT;
end
LM=sum(LocalScores(problem,:),2)/NumTrial1;
NLM=sum(NLocalScores(problem,:),2)/NumTrial1;
%h =heatmap(alpha_array,c_array,mean_precision);
%h.Title = ['Mean Precision, NoS: ', num2str(NumberOfSamples),' ', dataset,...
           %' Tau:',num2str(tau), ' Sigma:', num2str(sigma)];
%h.XLabel = 'NonLocality';
%h.YLabel = 'Teleportation';

Result(problem,1)=LM;
Result(problem,2)=NLM

end
save KFoldCVM Result LocalScores NLocalScores
