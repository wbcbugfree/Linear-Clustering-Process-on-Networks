function [comm,C,Q] = louvain(A)

A = relabel_graph(A);
A_ori = A;
comm = 1:size(A_ori,1);
uncomm = comm;
flag = false;
m = sum(A_ori, 'all')/2;%the sum of the weights of all links in the network

while flag == false
    flag = true;
    for i_c = 1:length(A) %the node that is about to move
        sigma_tot = sum(A);%sum of all edge weights for nodes within the community (including edges which link to other communities)
        sigma_tot = sigma_tot(:);
        sigma_tot(i_c) = 0;
        k_i = sum(A(:,i_c));%sum of the weights of the edges attached to nodes i 
        k_i_in = A(:,i_c);%the sum of the weights of the links between i and other nodes in the community that i is moving into
        k_i_in(i_c) = 0;
        diff_mod = k_i_in/m-(sigma_tot+k_i).*k_i/(2*m^2);
        [max_value, max_index]= max(diff_mod);
        if max_value>0
            flag = false;
            A(:,max_index) = A(:,max_index) + A(:,i_c);%moving nodes to maximum modularity and combine nodes
            A(max_index,:) = A(max_index,:) + A(i_c,:);
            A(:,i_c) = zeros(size(A,1),1);%removing node
            A(i_c,:) = zeros(1,size(A,2));
            tomove = uncomm(i_c);%assign communities to nodes
            dest = uncomm(max_index);
            comm(comm == tomove) = dest;
            uncomm(i_c) = 0;
        end
    end
    A = A((~uncomm==0),(~uncomm==0));
    uncomm = uncomm(~(uncomm==0));
end

uniComm = unique(comm);
for i=1:length(uniComm)
    comm(comm==uniComm(i)) = max(uniComm) +i;
end
uniComm = unique(comm);
for i=1:length(uniComm)
    comm(comm==uniComm(i)) = i;
end
Q = compute_modularity(A_ori,comm');
C = max(comm);
end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end

function A = relabel_graph(A)

G = graph(A);
edge = table2array(G.Edges);
new_edge = edge;
n = size(G.Nodes,1);

rand = randperm(n);
for i = 1:n
    new_edge(edge == i) = rand(i);
end

EdgeTable = table([new_edge(:,1) new_edge(:,2)]);
EdgeTable.Properties.VariableNames(1) = "EndNodes";
G = graph(EdgeTable);
A = adjacency(G);
end