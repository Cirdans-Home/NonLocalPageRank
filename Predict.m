function [predicted_edges,TP] = Predict(G,ind_deleted_edges,c,alpha,K,method,symmetry)
%Input:   G is the original graph
%         ind_deleted_edges indices of the edges to be deleted 
%         alpha NonLocal Parameter
%         c teleportation coefficient
%         K number of edges to be predicted
%         method is the method used to compute symilarity matrix
%         symmetry is True if a symmetric symilarity matrix is used
%Output:  Predicted Edges

n = numnodes(G);
e = ones(n,1);

H = G.rmedge(ind_deleted_edges);
AR = H.adjacency(); %Adjacency Matrix after removing edges
if ~isinf(alpha) && alpha > 0 
    % Nonlocal pagerank
    W = H.distances();
    W = 1./(W.^alpha);
    W(W == inf) = 0;
    D = (1./(W*e));
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    P = D*W;
elseif isinf(alpha)
    % Standard pagerank
    D = 1./(AR*e);
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    P = D*AR;
else
    error('alpha has to be 0 < \alpha <= \inf');
end
[X] = Similarity(P,method,c,symmetry);

if symmetry==true
    % We do not consider Diagonal
    T=triu(X,1);
    T(AR==1)=0;
    [~,J] = maxk(T(:),K);
    [i,j] = ind2sub([n n],J);
else
    % We do not consider Diagonal
    [~,J] = maxk(X-diag(diag(X)),K);
    [i,j] = ind2sub([n n],J);
end

predicted_edges=[i,j];

RE = G.Edges(ind_deleted_edges,:);
removed_edges=table2array(RE(:,1));
TP = sum(ismember(predicted_edges,removed_edges,'row')); 
end

