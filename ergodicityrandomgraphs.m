%% ERGODICITY FOR RANDOM GRAPHS
% ERGODICITY COEFFICIENTS FOR RANDOM GRAPHS computes the ergodicity
% coefficients for a the Small-World and Erdős-Rényi random graphs. Both
% the powerlaw and exponential transition can be selected for this task.
%
% The code makes use of functions from:
% Taylor, Alan, and Desmond J. Higham. "CONTEST: A controllable test matrix 
% toolbox for MATLAB." ACM Transactions on Mathematical Software (TOMS) 
% 35.4 (2009): 1-17.
% to generate the Small World and Erdős-Rényi random graphs.
%
% Code by:
% S. Cipolla - Università di Padova, Dipartimento di Matematica
% F. Durastante - Consiglio Nazionale delle Ricerche, Istituto per le
% Applicazioni del Calcolo "M. Picone"
% F. Tudisco - Gran Sasso Science Institute

clear; clc; close all;

rng(100);
typeoftransition = 'powerlaw'; %exponential
typeofgraph = 'erdos'; %smallworld
NumberOfGraphs = 100;   % Number of Random Graphs to be generated
NumberOfAlphas = 100;   % Number of Random \alphas to be generated
GraphSize = 500;        % Size of the graph, pay attention to the cubic cost
                        % of the tau1 computation!

for i=1:NumberOfGraphs
    
    switch typeofgraph
        case 'erdos'
            A = smallw(GraphSize,2,0.1);    % Using default options, see 
            % the CONTEST guide to use different settings
        case 'smallworld'
            A = erdrey(GraphSize);          % Using default options, see 
            % the CONTEST guide to use different settings
    end
    
    if issymmetric(A)
        G = graph(A,'omitselfloops');
    else
        G = digraph(A,'omitselfloops');
    end
    N = size(A,1);
    
    %% STANDARD PAGERANK
    e = ones(N,1);
    D = 1./(A*e);
    D(D == inf) = 0;
    Pinf = spdiags(D,0,N,N)*A;
    Ginf = 0.85*Pinf + (1-0.85)/N*(e*e.');
    erginf(i) = ergodicity(Ginf);
    
    %% NONLOCAL PAGERANK
    nodes = (1:N).';
    indexvalue = 1;
    distancematrix = G.distances();
    alphavalue = 5*rand(1,NumberOfAlphas);
    j = 1;
    for alpha = alphavalue
        switch typeoftransition
            case 'powerlaw'
                W = distancematrix;
                W = 1./(W.^alpha);
                W(W == inf) = 0;
                D = (1./(W*e));
                D(D == inf) = 0;
                D = spdiags(D,0,N,N);
                P = D*W;
            case 'exponential'
                W = distancematrix;
                W = exp(-alpha*W);
                W = sparse(W);
                W = W - spdiags(spdiags(W,0),0,N,N);
                D = (1./(W*e));
                D(D == inf) = 0;
                D = spdiags(D,0,N,N);
                P = D*W;
        end
        Galpha = 0.85*P + (1-0.85)/N*(e*e.');
        ergodicityvalue(i,j) = ergodicity(Galpha);
        j = j + 1;
    end
end

%% PLOT

figure(1)
semilogy(1:NumberOfGraphs,erginf,'o',1:NumberOfGraphs,ergodicityvalue,'kx')
xlabel('Test Graph');
ylabel('Ergodicity Coefficient');
legend({'Local PageRank','NonLocal PageRank'},'Location','southeast');
switch typeofgraph
    case 'erdos'
        title('Erdős-Rényi');
    case 'smallworld'
        title('Small World Network');
end
axis square

%% Auxiliary Functions from CONTEST
% These auxiliary functions come from
% Taylor, Alan, and Desmond J. Higham. "CONTEST: A controllable test matrix 
% toolbox for MATLAB." ACM Transactions on Mathematical Software (TOMS) 
% 35.4 (2009): 1-17.
% they are needed to generate the Small World and Erdős-Rényi random
% graphs.
function A = smallw(n,k,p)

%SMALLW     Generate adjacency matrix for a small world network.
%
%   Input   n: dimension of matrix (number of nodes in graph).
%           k: number of nearest-neighbours to connect. Defaults to 1.
%           p: probability of adding a shortcut in a given row. Defaults to
%           0.1.
%
%   Output  A: n by n symmetric matrix with the attribute sparse.
%
%   Description:    Shortcuts are added to a kth nearest neighbour ring
%                   network with n nodes by calling the utility function
%                   short.m.
%
%   Reference:  D.J. Watts, S. H. Strogatz,
%               Collective Dynamics of Small World Networks,
%               Nature 393 (1998), pp. 440-442.
%
%   Example:    A = smallw(100,1,0.2);

if nargin <= 2
    p = 0.1;
    if nargin == 1
        k = 2;
    end
end

twok = 2*k;

I = zeros(2*k*n,1);
J = zeros(2*k*n,1);
S = zeros(2*k*n,1);

for count = 1:n
    
    I( (count-1)*twok+1 : count*twok ) = count.*ones(twok,1);
    J( (count-1)*twok+1 : count*twok ) = mod([count:count+k-1 n-k+count-1:n+count-2],n)+1;
    S( (count-1)*twok+1 : count*twok ) = ones(twok,1);
    
end

A = sparse(I,J,S,n,n);
A = short(A,p);
end

function S = short(A,p)

%SHORT      Randomly add entries (shortcuts) to a matrix
%
%   Input   A: n by n adjacency matrix
%           p: probability that an entry is added to a given row
%
%   Output  S: n by n adjacency matrix with the attribute sparse.
%
%   Description:    A symmetric matrix of shortcuts is created which has
%                   an entry in each row with independent probability p.
%                   This is added to the matrix A.
%
%   Example: S = short(A,0.3);

n = length(A);

if nargin == 1
    p = log(n)/n;
end

Ihat = find(rand(n,1)<=p);
Jhat = ceil(n*rand(size(Ihat)));
Ehat = ones(size(Ihat));

self = find(Ihat==Jhat);
Ihat(self) = [];
Jhat(self) = [];
Ehat(self) = [];


[I,J,E] = find(A);

S = sparse([I;Ihat;Jhat],[J;Jhat;Ihat],[E;Ehat;Ehat],n,n);
end

function A = erdrey(n,m)

%ERDREY     Generate adjacency matrix for a G(n,m) type random graph.
%
%   Input   n: dimension of matrix (number of nodes in graph).
%           m: 2*m is the number of 1's in matrix (number of edges in graph).
%           Defaults to the smallest integer larger than n*log(n)/2.
%
%   Output  A: n by n symmetric matrix with the attribute sparse.
%
%
%   Description:    An undirected graph is chosen uniformly at random from
%                   the set of all symmetric graphs with n nodes and m
%                   edges.
%
%   Reference:  P. Erdos, A. Renyi,
%               On Random Graphs,
%               Publ. Math. Debrecen, 6 1959, pp. 290-297.
%
%   Example: A = erdrey(100,10);

if nargin == 1
    m = ceil(n*log(n)/2);
end

nonzeros = ceil(0.5*n*(n-1)*rand(m,1));
v = zeros(n,1);
for count = 1:n
    v(count) = count*(count-1)/2;
end

I = zeros(m,1);
J = zeros(m,1);
S = ones(m,1);

for count = 1:m
    i = min(find(v >= nonzeros(count)));
    j = nonzeros(count) - (i-1)*(i-2)/2;
    I(count) = i;
    J(count) = j;
end

A = sign(sparse([I;J],[J;I],[S;S],n,n));

while nnz(A) ~= 2*m
    
    difference = m-nnz(A)/2;
    Inew = zeros(difference,1);
    Jnew = zeros(difference,1);
    for count = 1:difference
        index = ceil(0.5*n*(n-1)*rand);
        Inew(count) = min(find(v>=index));
        Jnew(count) = index - (Inew(count)-1)*(Inew(count)-2)/2;
    end
    I = cat(1,I,Inew);
    J = cat(1,J,Jnew);
    S = ones(length(I),1);
    A = sign(sparse([I;J],[J;I],[S;S],n,n));
    
end
end