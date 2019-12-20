%%% SCRIPT TO CREATE THE MULTILAYER TENSOR OF METRO GRAPH 

% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute
clear all
clc
addpath('./tensor_toolbox-master');
A = dlmread('./London Tube Data/london_tube_edges.txt',' ',1,0);
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
