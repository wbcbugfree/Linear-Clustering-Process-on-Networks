function A = create_network( N, random_graph_model, p_ER, p_WS, K_WS, m0_BA, m_BA, connectivityFlag )
 switch random_graph_model
    case 'star'
            A = zeros( N );
            A( 1, : ) = 1;
            A( :, 1 ) = 1;
            A( 1, 1 ) = 0;
    case 'WS'
        G_graph = WattsStrogatz( N, K_WS, p_WS );
        A = G_graph.adjacency;
    case 'BA'
        A = BA( N, m0_BA,  m_BA );
    case 'ER_undirected'
        L = N*(N-1)/2;
        connectedComponents = 2;
        
        while connectedComponents>1 
            A = zeros(N);
            A(triu(ones(N), 1)>0) = (rand(L, 1) <= p_ER) ;
            A = (A + transpose(A));
            
            if connectivityFlag
                [connectedComponents, ~] = graphconncomp(sparse(A), 'Directed', false);
            else
                connectedComponents=1;
            end
        end
        
    case 'ER_directed'
        L = N*(N-1)/2;
        connectedComponents = 2;
        
        while connectedComponents>1
            A = zeros(N);
            A(triu(ones(N), 1)>0) = (rand(L, 1) <= p_ER) ;
            A(tril(ones(N), -1)>0) = (rand(L, 1) <= p_ER) ;
            
            if connectivityFlag
                [connectedComponents, ~] = graphconncomp(sparse(A), 'Directed', true);
            else
                connectedComponents=1;
            end
        end
end

if ~strcmp(random_graph_model, 'ER_directed') && norm( full( A - transpose( A )))>0
    error('Network not symmetric!')
end

end



function net = BA(  N, m0,  m  )
%Generates a network with N nodes and avg. node degree 2m following the
%Barabasi-Albert model, starting from a network with m0 nodes.
%
%This implementation is based on the model published in:
%
%   Barab?si,A.-L. and Albert,R. (1999) Emergence of Scaling in Random 
%   Networks. Science, 286, 509-512.
%
%INPUT
%   N -> Number of nodes in the graph (default = 1000).
%   m0 -> Initial number of nodes in the network (default = 5). This
%         parameters has to be > 0.
%   m -> Number of nodes with which a new node in the network connects to 
%        This parameter has to be <= to m0 (default = 3).
%
%OUPUT
%   net -> The adjacency matrix representation of the constructed 
%          Barb?si-Albert network. If the number of parameters exceeds the
%          maximum allowed, the returned matrix is empty.
%
% Copyright (C) Gregorio Alanis-Lobato, 2014

if m> m0
    error('m0 must be greater or equal than m')
end

%%% Barab?si-Albert network construction
net = zeros(N);
%Nodes 1 to m0 form a clique
net(1:m0, 1:m0) = 1;
net(logical(eye(N))) = 0;
for i = (m0+1):N
    %Add a new node to the network and connect to m existing nodes if they
    %indeed exist
    
    if i <= m
        net(i, 1:i) = 1; %If m is larger than the number of nodes in the system, connect to all of them
        net(i, i) = 0; %Self-connections not allowed
    else
        k = sum(net(1:(i-1), 1:(i-1)), 2); %Degree of the existing nodes
        k_total = sum(k);
        p = k./k_total; %Attraction probability of each node in the network
        
        r = rand(i-1, 1);
        
        if nnz(r < p) < m
            %Connect to the m minimum p
            [~, idx] = sort(r);
            net(i, idx(1:m)) = 1;
        else
            %Find with whom to connect to
            net(i, r < p) = 1;
        end     
        
    end
end
%Make sure the network is symmetric and has zero-diagonal
net = (net + transpose( net ))>0;
net(logical(eye(N))) = 0;
net = double(net);

[connectedComponents, ~] = graphconncomp(sparse(net), 'Directed', false);

if connectedComponents ~= 1
error('disconnected graph')
end

end

function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end

