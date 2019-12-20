function [D] = colored_distance_new(T)
% Input: T_ijk is a tensor of dimension n_nodes \times n_nodes \times n_kolors...
% such that T_ijk=1 if i~j with color k, T_ijk=0 otherwise 

% %%% Remark1: T must be such that T(:,:,k)=T(:,:,k)';
% %%% Remark2: we are assuming n_nodes >> n_kolors

% Output: matrix of ''colored distances''
N=size(T);
row_sums=zeros(N(1),N(3));
for k=1:N(3)
    M(N(1)*(k-1)+1:N(1)*k,N(1)*(k-1)+1:N(1)*k)=sparse(double(T(:,:,k)));
    row_sums(:,k)=sum(M(N(1)*(k-1)+1:N(1)*k,N(1)*(k-1)+1:N(1)*k),2);
end
row_sums=(row_sums~=0);
for i=1:N(3)-1
    for j=i+1:N(3)
        v=row_sums(:,i)&row_sums(:,j);
        v=double(v);
        M(N(1)*(i-1)+1:N(1)*i,N(1)*(j-1)+1:N(1)*j)=spdiags(v,0,N(1),N(1));
        M(N(1)*(j-1)+1:N(1)*j,N(1)*(i-1)+1:N(1)*i)=spdiags(v,0,N(1),N(1));
    end
end
A=graph(M);
DIST=A.distances();
D=Inf(N(1));
for i=1:N(3)
    for j=1:N(3)
        D=min(D,DIST(N(1)*(i-1)+1:N(1)*(i),N(1)*(j-1)+1:N(1)*j));
     end
end
end