%% SCRIPT FOR TABLES IN SECTION 5 

% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute
clear all
close all
clc
addpath('./tensor_toolbox-master');
% Script for  Tables 1-2

%%% The following line loads a Tensor (Tensor) and an adjacency Matrix (cMgraph) 
%%% representing 13 lines and 369 nodes
load('./london_tube_graphs.mat');
%%% Graph corresponding to 13 lines
G=graph(cMgraph);
%%% We want to consider just the subset of nodes which are on the lines 1-11
%%% These nodes coincide with the nodes which are in the file tube_usage_data.mat
load('London Tube Data/tube_usage_data.mat')
[usage_node_indexes,uni_perm]=sort(table2array(TABLE(:,1)),'ascend');
%% Ordering the table accordingly nodes numeration 
TABLE=TABLE(uni_perm,:);
%% Selecting just lines 1-11 and corresponding nodes
RTensor=Tensor(usage_node_indexes+1,usage_node_indexes+1,1:11);
RG=subgraph(G,usage_node_indexes+1);
%% Brief Recap on Problem Dimension
['Number of Nodes in lines 1--11 = ',num2str(size(usage_node_indexes,1))]
['Number of Edges in lines 1--11 = ',num2str(nnz(RG.adjacency)/2)]
CsY=zeros(12,10);
% %% Ground Truth Data Year 
for year=3:12
passengers=table2array(TABLE(:,year));
[gt,gtr]=sort(passengers,'descend');
ground_truth=usage_node_indexes(gtr);
%% Parameters
pr_alpha=0.85;
%% Classic PageRank
[score_pr] = test_pt(RG,distances(RG),Inf,pr_alpha);
[~,rank2]=sort(score_pr,'descend');
rank_pr=usage_node_indexes(rank2);
%% Different Ranking 
nonlocality_alpha=1.3;
[D] = colored_distance_new(RTensor);

% Metro Distance used
%[score_md] = test_pt(RG,D,nonlocality_alpha,pr_alpha);
[score_md] = test_pt(RG,D,1.7,pr_alpha);
[~,rank1_md]=sort(score_md,'descend');
rank_md=usage_node_indexes(rank1_md);

% Shortest Path Distance
[score_spd] = test_pt(RG,distances(RG),1.7,pr_alpha);
[~,rank1_spd]=sort(score_spd,'descend');
rank_spd=usage_node_indexes(rank1_spd);

%% Computing Cumlative Sums
cs_gt=cumsum(passengers(gtr));
cs_pr=cumsum(passengers(rank2));
cs_md=cumsum(passengers(rank1_md));
cs_spd=cumsum(passengers(rank1_spd));

% % CsY(:,year-2)=[cs_gt(5);cs_gt(10);cs_gt(15);cs_gt(20);...
% %     cs_pr(5);cs_pr(10);cs_pr(15);cs_pr(20);...
% %     cs_spd(5);cs_spd(10);cs_spd(15);cs_spd(20);...
% %     cs_md(5);cs_md(10);cs_md(15);cs_md(20)];
% 
CsY(:,year-2)=[cs_gt(5);cs_gt(15);cs_gt(45);...
    cs_pr(5);cs_pr(15);cs_pr(45);...
    cs_spd(5);cs_spd(15);cs_spd(45);...
    cs_md(5);cs_md(15);cs_md(45)];
end 
save('Millions_passengers_results.mat','CsY');
T=array2table(CsY);
table2latex(T,'tenyear')
%ltable=latex(vpa(sym(CsY),5))

%%% Tables

%% Ground Truth Data Year 
for year=3:3
passengers=table2array(TABLE(:,year));
[gt,gtr]=sort(passengers,'descend');
ground_truth=usage_node_indexes(gtr);
%% Parameters
pr_alpha=0.85;
%% Classic PageRank
[score_pr] = test_pt(RG,distances(RG),Inf,pr_alpha);
[~,rank2]=sort(score_pr,'descend');
rank_pr=usage_node_indexes(rank2);
nonlocality_alpha=1.7;
[D] = colored_distance_new(RTensor);

% Metro Distance used
[score_md] = test_pt(RG,D,nonlocality_alpha,pr_alpha);
[~,rank1_md]=sort(score_md,'descend');
rank_md=usage_node_indexes(rank1_md);

% Shortest Path Distance
[score_spd] = test_pt(RG,distances(RG),nonlocality_alpha,pr_alpha);
[~,rank1_spd]=sort(score_spd,'descend');
rank_spd=usage_node_indexes(rank1_spd);

%Ground_Truth=[TABLE.(2)(gtr(1:10)),TABLE.(3)(gtr(1:10))];
% Local_PageRank=TABLE(rank2(1:10),2:3);
% NonLocal_PageRankSP=TABLE(rank1_spd(1:10),2:3);
% NonLocal_PageRankMD=TABLE(rank1_md(1:10),2:3);
% T = table(Ground_Truth,Local_PageRank,NonLocal_PageRankSP,NonLocal_PageRankMD);
T = table(TABLE.(2)(gtr(1:10)),TABLE.(3)(gtr(1:10)),TABLE.(2)(rank2(1:10)),TABLE.(3)(rank2(1:10)),...
    TABLE.(2)(rank1_spd(1:10)),TABLE.(3)(rank1_spd(1:10)),TABLE.(2)(rank1_md(1:10)),TABLE.(3)(rank1_md(1:10)))
%T = table(TABLE(gtr(1:10),2:3),TABLE(rank2(1:10),2:3),TABLE(rank1_spd(1:10),2:3),TABLE(rank1_md(1:10),2:3));
end
save('Millions_passengers_results1.mat','T');
save('vettoriperpython.mat','ground_truth', 'rank_pr','rank_md', 'rank_spd'); 
table2latex(T,'topten')
%ltable=latex(vpa(sym(CsY),5))


