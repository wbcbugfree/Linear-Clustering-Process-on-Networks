function [Comm,C,Q] = Leiden(A)
A_0 = A;                                % Store the original adjacency matrix
N = size(A,1);                          % Number of nodes
P = 1:N;                                % Each node is a community
done = 0;                               % Indicator for stopping the code
it = 1;                                 % Initialise the iteration counter

while ~done
    P = MoveNodesFast(A,P);
    P_refined = RefinePartition(A,P);
    P = simplify_C(P_refined);                  % Define partition from 1 to c
    C_st{1,it} = P;                     % Store the current partition before aggregaton
    done = (length(unique(P)) == N);    % Check if there were no changes
    if(~done)                           % In case there were changes
        [A,P,N] = AggregateGraph(A,P);  % Aggregate the graph in c nodes
        it = it + 1;                    % Increment the iteration counter
    end
end
Comm = compute_partition(C_st,size(A_0,1));% Compute the final partition
C = length(unique(Comm));                  % Compute the number of clusters c
Q = compute_modularity(A_0,Comm);          % Compute modularity
end

function P = MoveNodesFast(A,P)
N = size(A,1);          
queue = randperm(N);                                        % Visit nodes at random
L_2 = sum(sum(A));                                          % Sum of all link weights
while ~isempty(queue)
    node_v = queue(1);
    queue(1) = [];
    ind_neig = intersect(find(P~=P(node_v)),find(A(:,node_v)));        % Find neighbours of node v
    Delta_m = zeros(length(ind_neig),1);                    % Initialise the delta_m vector of node v
    for counter_neig = 1 : length(ind_neig)                 % For each nieghbour j of node v
        ind_comm = find(P == P(ind_neig(counter_neig)));    % Determine all nodes from the community of node j
        Sum_in = sum(sum(A(ind_comm,ind_comm)));            % Sum of links within the community C
        Sum_tot = 2*sum(sum(A(ind_comm,:)));                % Sum of the weights of links incident to nodes in community C
        Sum_C = sum(A(ind_comm,node_v));            % Sum of the weights of links from node v to nodes in community C
        d_i = sum(A(:,node_v));                     % Sum of the link weights adjacent to node v
        Delta_m(counter_neig) = ((Sum_in + 2*Sum_C)/L_2 - ((Sum_tot + d_i)/L_2)^2) - (Sum_in/L_2 - (Sum_tot/L_2)^2 - (d_i/L_2)^2);   % Compute the modularity gain
    end
    [ind_val,ind_pos] = max(Delta_m);                       % Determine the best community for node i
    if ind_val > 0                                      
        P(node_v) = P(ind_neig(ind_pos));        % Update the cluster membership of each node
        for i = 1:length(ind_neig)
            if P(ind_neig(i))~=P(node_v)&&~ismember(ind_neig(i),queue)
                queue = [queue ind_neig(i)];
            end
        end
        A(:,ind_neig(ind_pos)) = A(:,ind_neig(ind_pos)) + A(:,node_v);%moving nodes to maximum modularity and combine nodes
        A(ind_neig(ind_pos),:) = A(ind_neig(ind_pos),:) + A(node_v,:);
        A(:,node_v) = zeros(size(A,1),1);%removing node
        A(node_v,:) = zeros(1,size(A,2));
    end
end
end

function P_refined = RefinePartition(A,P)
N = size(A,1);
P_refined = 1:N;
C = unique(P);
for i = 1:length(C)
    S = find(P == C(i));
    P_refined = MergeNodesSubset(A,P_refined,S,C(i));
end

end

function P = MergeNodesSubset(A,P,S,C)
R = [];
gamma = 1/7;
for i = 1:length(S)
    S_v = setdiff(S,S(i));
    recur_v = sum(A(:,S(i)));
    recur_S = length(S);
    E_v_Sv = length(intersect(find(A(:,S(i))), S_v));
    if E_v_Sv >= gamma * recur_v * (recur_S - recur_v)
        R = [R S(i)];
    else
        P(S(i)) = C;
    end
end

theta = 1;
for i = 1:length(R)
    T = [];
    if R(i)==P(R(i))
        for j = 1:length(S)
            S_C = setdiff(S,S(j));
            recur_C = sum(A(:,S(j)));
            recur_S = length(S);
            E_C_SC = length(intersect(find(A(:,S(j))), S_C));
            if E_C_SC >= gamma * recur_C * (recur_S - recur_C)
                T = [T S(j)];
            end
        end
        L_2 = sum(sum(A)); 
        Delta_m = zeros(length(T),1);                    % Initialise the delta_m vector of node v
        for counter_T = 1 : length(T)                 % For each nieghbour j of node v
            ind_comm = find(P == P(T(counter_T)));    % Determine all nodes from the community of node j
            Sum_in = sum(sum(A(ind_comm,ind_comm)));            % Sum of links within the community C
            Sum_tot = 2*sum(sum(A(ind_comm,:)));                % Sum of the weights of links incident to nodes in community C
            Sum_C = sum(A(ind_comm,R(i)));            % Sum of the weights of links from node v to nodes in community C
            d_i = sum(A(:,R(i)));                     % Sum of the link weights adjacent to node v
            Delta_m(counter_T) = ((Sum_in + 2*Sum_C)/L_2 - ((Sum_tot + d_i)/L_2)^2) - ...
                (Sum_in/L_2 - (Sum_tot/L_2)^2 - (d_i/L_2)^2);   % Compute the modularity gain
        end
        if Delta_m > 0
            Pr_C = exp((1/theta) * Delta_m);
            Pr_C = Pr_C/sum(Pr_C);
            Random = randsrc(1,1,[1:length(Pr_C); Pr_C']);
            P(R(i)) = T(Random);
        else
            P(R(i)) = C;
        end
    else
        P(R(i)) = C;
    end
end
end
    
function Q = compute_modularity(A,C) % Compute modularity of the given partition
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*C - C'*ones(1,N)) == 0)));
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
