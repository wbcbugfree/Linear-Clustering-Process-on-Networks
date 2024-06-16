function [C,c,Q] = Louvain_Algorithm_v2(A)
A_0 = A;                                % Store the original adjacency matrix
N = size(A,1);                          % Number of nodes
P = 1:N;                                % Each node is a community
Q = compute_modularity(A_0,P);          % Compute modularity
done = 0;                               % Indicator for stopping the code
it = 1;                                 % Initialise the iteration counter

while(~done)                            % Until modulairty m cannot be improved further
    P = MoveNodes(A,P,Q);               % Move nodes to the best community
    P = simplify_C(P);                  % Define partition from 1 to c
    C_st{1,it} = P;                     % Store the current partition before aggregaton
    Q = compute_modularity(A,P);        % Compute modularity
    done = (length(unique(P)) == N);    % Check if there were no changes
    if(~done)                           % In case there were changes
        [A,P,N] = AggregateGraph(A,P);  % Aggregate the graph in c nodes
        Q = compute_modularity(A,P);    % Compute modularity
        it = it + 1;                    % Increment the iteration counter
    end
end
C = compute_partition(C_st,size(A_0,1));% Compute the final partition
c = length(unique(C));                  % Compute the number of clusters c
Q = compute_modularity(A_0,C);          % Compute modularity
end

function P = MoveNodes(A,P,m_old)
N = size(A,1);          
m = m_old;
Deg = A*ones(N,1);
Q_p = 1/nnz(A).*(A - (1./nnz(A).*Deg*Deg'));
while(m >= m_old)
    nodes = randperm(N);                                        % Visit nodes at random
    L_2 = sum(sum(A));                                          % Sum of all link weights
    for counter = 1 : N                                         % Go over each node in the graph
        ind_neig = intersect(find(P~=P(nodes(counter))),find(A(:,nodes(counter))));        % Find neighbours of node i
        Delta_m = zeros(length(ind_neig),1);                    % Initialise the delta_m vector of node i
        for counter_neig = 1 : length(ind_neig)                 % For each nieghbour j of node i
            ind_comm = find(P == P(ind_neig(counter_neig)));    % Determine all nodes from the community of node j
            Sum_in = sum(sum(A(ind_comm,ind_comm)));            % Sum of links within the community C
            Sum_tot = 2*sum(sum(A(ind_comm,:)));                % Sum of the weights of links incident to nodes in community C
            Sum_C = sum(A(ind_comm,nodes(counter)));            % Sum of the weights of links from node i to nodes in community C
            d_i = sum(A(:,nodes(counter)));                     % Sum of the link weights adjacent to node i
            Delta_m(counter_neig) = ((Sum_in + 2*Sum_C)/L_2 - ((Sum_tot + d_i)/L_2)^2) - (Sum_in/L_2 - (Sum_tot/L_2)^2 - (d_i/L_2)^2);   % Compute the modularity gain
        end
        [ind_val,ind_pos] = max(Delta_m);                       % Determine the best community for node i
        if(ind_val > 0)                                         % If the modularity can be improved?
            ind_sub_1 = find(P == P(nodes(counter))); ind_sub_1(ind_sub_1~=nodes(counter));
            P(nodes(counter)) = P(ind_neig(ind_pos));        % Update the cluster membership of each node
            ind_sub_2 = find(P == P(ind_neig(ind_pos))); ind_sub_2(ind_sub_2~=nodes(counter));
            m = m - 2*sum(Q_p(nodes(counter),ind_sub_1)) + 2*sum(Q_p(nodes(counter),ind_sub_2));
        end
    end
    if(m <= m_old)                                              % If the modularity is not improved
        return                                                  % Break the while loop
    end
    m_old = m;                                                  % Update the old modularity
end
end
    
function [A_new,P_new,N_new] = AggregateGraph(A,P) % Aggregate the graph
N_new = length(P);
A_new = zeros(N_new);
P_new = 1 : N_new;
for counter_i = 1 : N_new
    for counter_j = 1 : N_new
        if(counter_i == counter_j)
            ind_comm_i = find(P == counter_i);
            A_new(counter_i,counter_i) = 0.5*sum(sum(A(ind_comm_i,ind_comm_i)));
        else
            ind_comm_i = (P == counter_i);
            ind_comm_j = (P == counter_j);
            A_new(counter_i,counter_j) = sum(sum(A(ind_comm_i,ind_comm_j)));
        end
    end
end
end

function Q = compute_modularity(A,C) % Compute modularity of the given partition
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*C - C'*ones(1,N)) == 0)));
end

function  P = simplify_C(P) % Define clusters from 1 to c
P_un = unique(P);
for counter = 1 : length(P_un)
    ind = (P == P_un(counter));
    P(ind) = counter;
end
end

function C = compute_partition(C_st,N) % Restore the partition from C_st
it = length(C_st);
C = zeros(1,N);
for counter_1 = 1 : N
    C(counter_1) = C_st{1,1}(counter_1);
    for counter_2 = 2 : it
        C(counter_1) = C_st{1,counter_2}(C(counter_1));
    end
end
end
