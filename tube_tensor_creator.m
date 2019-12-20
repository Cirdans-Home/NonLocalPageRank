%type london_tube_multiplex.edges
%%% This script creates the mathematical objects needed to analize London
%%% tube. The output includes 11 metro lines + 2 train lines and all 
%%% corresponding nodes.
clear all
clc
addpath('./tensor_toolbox-master');
%fileID = fopen('./datasets/London_Multiplex_Transport/Dataset/london_tube_multiplex.edges','r');
%formatSpec = '%d';
%A = textscan(fileID,formatSpec)
A = dlmread('./datasets/London_Multiplex_Transport/Dataset/london_tube_multiplex_new.edges',' ',1,0);
%A=int64(A);
%fclose(fileID);
[n,m]=size(A);
maxv=max(A);
modes=maxv(1);
nodes=max(maxv(2:3))+1;
Tensor = zeros(nodes,nodes,modes);
% In Each unfolding the graph of line of the metro
for j=1:n
    if A(j,2)+1~=A(j,3)+1
    Tensor(A(j,2)+1,A(j,3)+1,A(j,1))=1;
    Tensor(A(j,3)+1,A(j,2)+1,A(j,1))=1;
    end
end
% Graph 
% cMgraph=zeros(nodes,nodes);
% for j=1:modes
%     cMgraph=cMgraph+Tensor(:,:,j);
% end

cMgraph=sum(Tensor,3);
cMgraph=sparse(sign(cMgraph));

save 'london_tube_graphs.mat' Tensor cMgraph