%% SCRIPT FOR FIGURE 9

% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute
clear all
close all
clc
% Script for Figure 9
addpath('./tensor_toolbox-master');
%%% This load a Tensor (Tensor) and an adjacency Matrix (cMgraph) 
%%% representing 13 lines and 369 nodes
load('./london_tube_graphs.mat');
%%% Graph corresponding to 13 lines
G=graph(cMgraph);
%%% We want to consider just the subset of nodes which are on the lines 1-11
%%% These nodes coincide with the nodes which are in the file tube_usage_data.mat
load('./datasets/London_Multiplex_Transport/Tube_usage_data/tube_usage_data.mat')
[usage_node_indexes,uni_perm]=sort(table2array(TABLE(:,1)),'ascend');
%% Ordering the table accordingly nodes numeration 
TABLE=TABLE(uni_perm,:);
%% Selecting just lines 1-11 and corresponding nodes
RTensor=Tensor(usage_node_indexes+1,usage_node_indexes+1,1:11);
RG=subgraph(G,usage_node_indexes+1);
%% Brief Recap on Problem Dimension
['Number of Nodes in lines 1--11 = ',num2str(size(usage_node_indexes,1))]
['Number of Edges in lines 1--11 = ',num2str(nnz(RG.adjacency)/2)]
%% Ground Truth Data Year 2017
passengers=table2array(TABLE(:,3));
[gt,gtr]=sort(passengers,'descend');
ground_truth=usage_node_indexes(gtr);
%% PageRank+Parameters
pr_alpha=0.85;
[score_pr] = test_pt(RG,distances(RG),Inf,pr_alpha);
[~,rank2]=sort(score_pr,'descend');
rank_pr=usage_node_indexes(rank2);
[pr_gt,~]=isim_new(rank_pr,ground_truth);
cs_pr=cumsum(passengers(rank2));
%% Testing Ranking For Different Alphas
alpha_array=0.1:0.2:5;
[D] = colored_distance_new(RTensor);
for i=1:size(alpha_array,2)
% Metro Distanceused
[score_md(:,i)] = test_pt(RG,D,alpha_array(i),pr_alpha);
[~,rank1_md]=sort(score_md(:,i),'descend');
rank_md(:,i)=usage_node_indexes(rank1_md);
% Shortest Path Distance
[score_spd(:,i)] = test_pt(RG,distances(RG),alpha_array(i),pr_alpha);
[~,rank1_spd]=sort(score_spd(:,i),'descend');
rank_spd(:,i)=usage_node_indexes(rank1_spd);
% Computing Isim
[md_gt(:,i),~]=isim_new(rank_md(:,i),ground_truth);
[spd_gt(:,i),~]=isim_new(rank_spd(:,i),ground_truth);
% Computing Ratios
% Metro Distance
isim_ratios_md_gt(:,i)=md_gt(:,i)./pr_gt;
cs_md(:,i)=cumsum(passengers(rank1_md));
cs_ratios_md(:,i)=cs_md(:,i)./cs_pr;
% Shortest Path
isim_ratios_spd_gt(:,i)=spd_gt(:,i)./pr_gt;
cs_spd(:,i)=cumsum(passengers(rank1_spd));
cs_ratios_spd(:,i)=cs_spd(:,i)./cs_pr;
end 
save 'tube_results.mat'

%%
% figure(1)
% plot(alpha_array,isim_ratios_spd_gt(5,:),'k-o','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,isim_ratios_spd_gt(10,:),'k-*','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,isim_ratios_spd_gt(15,:),'k-x','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,isim_ratios_md_gt(5,:),'r-o','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,isim_ratios_md_gt(10,:),'r-*','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,isim_ratios_md_gt(15,:),'r-x','LineWidth',2, 'Markersize',10)
% legend('Spd-Gt 5','Spd-Gt 10','Spd-Gt 15','Md-Gt 5','Md-Gt 10','Md-Gt 15','Location', 'Best');
% title('Isim Ratios')
% figure(2)
% plot(alpha_array,cs_ratios_spd(5,:),'k-o','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,cs_ratios_spd(10,:),'k-*','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,cs_ratios_spd(15,:),'k-x','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,cs_ratios_md(5,:),'r-o','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,cs_ratios_md(10,:),'r-*','LineWidth',2, 'Markersize',10)
% hold on
% plot(alpha_array,cs_ratios_md(15,:),'r-x','LineWidth',2, 'Markersize',10)
% legend('Spd-Gt 5','Spd-Gt 10','Spd-Gt 15','Md-Gt 5','Md-Gt 10','Md-Gt 15','Location', 'Best');
% title('Cumsum Ratios')

%% Figures 
%%
id = 3:13;
figure

subplot(131)
plot(alpha_array(id),isim_ratios_spd_gt(15,id),'-','LineWidth',2, 'Markersize',10)
hold on
plot(alpha_array(id),isim_ratios_md_gt(15,id),'-','LineWidth',4, 'Markersize',10)
xline(alpha_array(9), '--','linewidth',1)
set(gca,'fontsize',12)
legend('SP distance','Multilayer distance','Location', 'northeast');
ylabel('Top 15 nodes ISIM ratio vs local PageRank','fontsize',15);
xlabel('Nonlocality parameter \alpha', 'fontsize',15);
xlim([0.5 2.5])
axis square

subplot(132)
plot(alpha_array(id),cs_ratios_spd(15,id),'-','LineWidth',2, 'Markersize',10)
hold on
plot(alpha_array(id),cs_ratios_md(15,id),'-','LineWidth',4, 'Markersize',10)
xline(alpha_array(9), '--','linewidth',1)
set(gca,'fontsize',12)
legend('SP distance','Multilayer distance','Location', 'southeast');
ylabel('Top 15 nodes cumsum ratio vs local PageRank','fontsize',15);
xlabel('Nonlocality parameter \alpha', 'fontsize',15);
xlim([0.5 2.5])
axis square


subplot(133)
range = 1:25;
plot(spd_gt(range,9),'-','LineWidth',2, 'Markersize',10)
hold on
plot(md_gt(range,9),'-','LineWidth',4, 'Markersize',10)
plot(pr_gt(range),'-','LineWidth',2, 'Markersize',10)
set(gca,'fontsize',12)
legend('SP distance','Multilayer distance','Local PageRank','Location', 'southeast');
ylabel('ISIM comparison for \alpha=1.7','fontsize',15);
xlabel('# top ranked nodes', 'fontsize',15);
axis square
set(gcf,'color','w');
tightfig
savefigures = true;
keyboard
if savefigures   
     filename = './result_figures/isim_comparison';
     savefig(filename);
     set(gcf, 'PaperPositionMode', 'auto');
     print(filename,'-depsc2', '-r600'); 
end
