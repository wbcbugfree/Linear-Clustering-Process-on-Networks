function [A,C,c,m] = LFR_benchmark_v9(N,d_av,c,gamma,beta,mu)
% N     - number of nodes
% d_av  - average degree
% c     - number of communities
% gamma - power law exponent of the degree distribution
% beta  - power law exponent of the community size distribution
% mu    - average ratio of inter-community links per node
clc
close all
dmin = 5;                                                   % Define the minimum degree value
dmax = 100;                                                 % Define the maximum degree value
Deg = powersample_Deg(gamma, N, dmin, dmax, d_av*N);        % Compose the degree vector
A = graphFromDeg(Deg);                                      % Compose the graph from the degree vector
N_comm = powersample_Com(beta, c, dmin + 1, dmax - 1, N);   % Compose the community size vector
C = assign_community(A,N_comm);                             % Assign a community to each node
A = rewire_graph(A,C,N_comm,mu);                            % Rewire the graph to obey mu rule
m = compute_modularity(A,C');                               % Compute the modularity

end

%% Used functions
function x = powersample_Com(beta, N, xmin, xmax, xsum) % Community size vector distribution obeying power law
% Generate power-law distributed variables
% gamma: power-law exponent
% N: graph size
% c: number of samples
% xmin: minimum value
% xmax: maximum value

% Generate a uniform sample in [0,1]
u = rand(N, 1);
% Calculate the corresponding power-law distributed variable
x = (xmax^(-beta+1) - xmin^(-beta+1)) * u + xmin^(-beta+1);
x = x.^(1/(-beta+1));
% Adjust values that are out of range
x(x < xmin) = xmin;
x(x > xmax) = xmax;

x = x/sum(x)*xsum;
x = ceil(x);
while(sum(x)~= xsum)
    ind = ceil(rand*N);
    if(sum(x) > xsum && x(ind) > xmin)
        x(ind) = x(ind) - 1;
    elseif(sum(x) < xsum && x(ind) < xmax)
        x(ind) = x(ind) + 1;
    end
end
end

function d = powersample_Deg(gamma, N, dmin, dmax, L2) % Degree vector distribution obeying power law
% Generate power-law distributed variables
% gamma: power-law exponent
% N: number of samples (graph size)
% dmin: minimum degree value
% dmax: maximum degree value
% L2: sum of values (number of links)

% Generate a uniform sample in [0,1]
u = rand(N, 1);
% Calculate the corresponding power-law distributed variable
d = (dmax^(-gamma+1) - dmin^(-gamma+1)) * u + dmin^(-gamma+1);
d = d.^(1/(-gamma+1));
% Adjust values that are out of range
d(d < dmin) = dmin;
d(d > dmax) = dmax;
d = d/sum(d)*L2;
d = ceil(d);
if(rem(sum(d),2))
    ind = find(rem(d,2) == 1);
    d(ind(1)) = d(ind(1)) + 1;
end
end

function A = graphFromDeg(Deg) % Construct the graph from the degree sequence
A = zeros(length(Deg));
while sum(Deg) > 0  % while there are still stubs to connect
    % order stubs by decreasing number of degrees left
    [sorted,I] = sort(-Deg);
    n1 = I(1);
    for x=1:-sorted(1)
        n2 = I(x+1);
        A(n1,n2)= A(n1,n2)+1;
        A(n2,n1)= A(n2,n1)+1;
        Deg(n1) = Deg(n1)-1;
        Deg(n2) = Deg(n2)-1;
    end
end
end

function C = assign_community(A,N_comm) % Assign community to each node
c = length(N_comm);                 
N = size(A,1);
Pool = 1:N;
C = zeros(N,1);
while(~isempty(find(C == 0, 1)))
    ind_node = ceil(rand*length(Pool));
    ind_neig = A(ind_node,:) == 1;
    flag = 1;
    while(flag)
        ind_comm = ceil(rand*c);
        int_deg = length(find(C(ind_neig) == ind_comm));
        if(N_comm(ind_comm) > int_deg && length(find(C == ind_comm)) < N_comm(ind_comm))
            C(ind_node) = ind_comm;
            flag = 0;
        end
    end
end
end

function A = rewire_graph(A,C,N_comm,mu) % Rewire the graph, to obey the mixing constraint
N = size(A,1);                                          % Store the graph size
Deg = A*ones(N,1);                                      % Compute the degree vector
IN = compute_IN(A,Deg,C);                               % Compute mu per each node
B = compute_B(A,C,N_comm,IN);                           % Compute the modified incidence matrix
cnt_0 = 0;                                              % Total number of reqiring counter
cnt = 1;                                                % Succesive failed rewiring counter
while((mean(IN) > 1.1*mu || mean(IN) < 0.90*mu || var(IN) > 0.009) && (cnt < 10) && (cnt_0 < nnz(A)/4))
    cnt_0 = cnt_0 + 1;                                  % Update the total number of rewiring counter
    if(mean(IN) > 1.1*mu)                               % If the average mu is too high 
        [n_1,n_2,n_3,n_4] = decrease_mu_v2(A,B,mu);     % Decrease mu by rewiring a pair of links
        if(n_4 > 0)                                     % If the rewiring was succesful
            A = rewireG(A,n_1,n_2,n_3,n_4);             % Rewire the graph
            IN = compute_IN(A,Deg,C);                   % Recompute the mu vector
            B = compute_B(A,C,N_comm,IN);               % Recompute the B matrix
            cnt = 1;                                    % Reset the succesive failed rewiring counter
        else
            cnt = cnt + 1;                              % Else, increment the succesive failed rewiring counter
        end
    elseif(mean(IN) < 0.9*mu)                           % If the average mu is too low
        [n_1,n_2,n_3,n_4] = increase_mu(A,B,mu);        % Increase mu by rewiring a pair of links
        if(n_4 > 0)                                     % If the rewiring was succesful
            A = rewireG(A,n_1,n_2,n_3,n_4);             % Rewire the graph
            IN = compute_IN(A,Deg,C);                   % Recompute the mu vector
            B = compute_B(A,C,N_comm,IN);               % Recompute the B matrix
            cnt = 1;                                    % Reset the succesive failed rewiring counter
        else
            cnt = cnt + 1;                              % Else, increment the succesive failed rewiring counter                     
        end
    end
    if(max(IN) > 1.2*mu)                                % If the maximum mu is too high
        [n_1,n_2,n_3,n_4] = decrease_mu_v2(A,B,mu);     % Decrease mu by rewiring a pair of links
        if(n_4 > 0)                                     % If the rewiring was succesful
            A = rewireG(A,n_1,n_2,n_3,n_4);             % Rewire the graph
            IN = compute_IN(A,Deg,C);                   % Recompute the mu vector
            B = compute_B(A,C,N_comm,IN);               % Recompute the B matrix
            cnt = 1;                                    % Reset the succesive failed rewiring counter
        else
            cnt = cnt + 1;                              % Else, increment the succesive failed rewiring counter
        end
    end
    if(min(IN) < 0.8*mu)                                % If the minimum mu is too low
        [n_1,n_2,n_3,n_4] = increase_mu(A,B,mu);        % Increase mu by rewiring a pair of links
        if(n_4 > 0)                                     % If the rewiring was succesful
            A = rewireG(A,n_1,n_2,n_3,n_4);             % Rewire the graph
            IN = compute_IN(A,Deg,C);                   % Recompute the mu vector
            B = compute_B(A,C,N_comm,IN);               % Recompute the B matrix
            cnt = 1;                                    % Reset the succesive failed rewiring counter
        else
            cnt = cnt + 1;                              % Else, increment the succesive failed rewiring counter
        end
    end
end
clc
%disp('Graph generated!')
%disp('Total number of iterations:')
%disp(cnt_0)
end

function [n_1,n_2,n_3,n_4] = decrease_mu_v2(A,B,mu) % Increase the ratio of inter-community links in the graph
% This function identifies best two pair of connected nodes to perform
% rewiring with. Since the aim is to decrease mu, we want to replace two
% inter-community links with two intra-community links.
ind = intersect(find(B(3,:) > mu),find(B(4,:) > mu));               % Identify only those nodes with IN value larger than 1.1*mu
B = B(:,ind);
ind = find(B(5,:) - B(6,:) ~= 0);                                   % Identify only inter-community links
B = B(:,ind);
B = B(:,randperm(length(ind)));
cnt = 1;
while(cnt < length(ind))
    n_1 = B(1,cnt);
    n_2 = B(2,cnt);
    Com_1 = B(5,cnt);
    Com_2 = B(6,cnt);
    ind_temp = union(intersect(find(B(5,:) == Com_1),find(B(6,:) == Com_2)),...
                     intersect(find(B(5,:) == Com_2),find(B(6,:) == Com_1)));
    B_temp = B(:,ind_temp);
    for counter = 1 : length(ind_temp)
        if(B_temp(5,counter) == Com_1)
            n_3 = B_temp(1,counter);
            n_4 = B_temp(2,counter);
        else
            n_4 = B_temp(1,counter);
            n_3 = B_temp(2,counter);
        end
        if(A(n_1,n_3)==0 && A(n_1,n_4)==0 && A(n_2,n_3)==0 && A(n_2,n_4)==0)
            return
        end
    end
    cnt = cnt + 1;
end
n_1 = 0; n_2 = 0; n_3 = 0; n_4 = 0;
end

function [n_1,n_2,n_3,n_4] = increase_mu(A,B,mu) % Increase the ratio of inter-community links in the graph
% This function identifies best two pair of connected nodes to perform
% rewiring with. Since the aim is to increase mu, we want to replace two
% inter-community links with two intra-community links. Therefore, we
% identify two pairs of connected nodes from the same community.
ind = intersect(find(B(3,:) < mu),find(B(4,:) < mu));       % Identify only those nodes with IN value larger than 1.1*mu
B = B(:,ind);
%% Identify the first pair of nodes
ind = find(B(5,:) - B(6,:) == 0);                               % Identify only intra-community links
B = B(:,ind);
[~,pos_IN] = sort(B(3,:) + B(4,:),'ascend');                    % Sort the intra-community links based on their IN value, in descending order
pos = ceil(rand*0.4*length(pos_IN));
if(isempty(ind))
    n_1 = 0; n_2 = 0; n_3 = 0; n_4 = 0;
    return
end
n_1 = B(1,(pos_IN(pos)));                                       % Choose the link with the highest IN value
n_2 = B(2,(pos_IN(pos)));                                       % Choose the second link
pos_12 = (pos_IN(pos));                                         % Store the position of the link (n_1 ~ n_2)
Com_12 = B(5,pos_12);
%% Exclude neighbouring nodes of node 1 and node 2 from consideration
ind_neig_1 = find(A(:,n_1) == 1);
ind_neig_2 = find(A(:,n_2) == 1);
% Look for another community
ind = find(B(5,:) ~= Com_12);
B = B(:,ind);
for counter_1 = 1 : length(ind_neig_1) % Remove links adjacent to neighbouring nodes of node 1
    ind = intersect(find(B(1,:) ~= ind_neig_1(counter_1)),find(B(2,:) ~= ind_neig_1(counter_1)));
    B = B(:,ind);
end
for counter_2 = 1 : length(ind_neig_2) % Remove links adjacent to neighbouring nodes of node 2
    ind = intersect(find(B(1,:) ~= ind_neig_2(counter_2)),find(B(2,:) ~= ind_neig_2(counter_2)));
    B = B(:,ind);
end
if(isempty(ind))
    n_3 = 0; n_4 = 0;
    return
end
%% Identify the second pair of nodes
[~,pos_IN] = sort(B(3,:) + B(4,:),'ascend');
n_3 = B(1,pos_IN(1));
n_4 = B(2,pos_IN(1));
end

function IN = compute_IN(A,Deg,C) % COmpute the ratio of inter-community links per node
N = size(A,1);
IN = zeros(N,1);
for counter = 1 : N
    IN(counter) = length(intersect(find(A(counter,:) == 1),find(C ~= C(counter))))/Deg(counter);
end
end

function A = rewireG(A,ind_node_1,ind_node_2,ind_node_3,ind_node_4) % Perform rewiring
A(ind_node_1,ind_node_2) = 0;
A(ind_node_2,ind_node_1) = 0;

A(ind_node_3,ind_node_4) = 0;
A(ind_node_4,ind_node_3) = 0;

A(ind_node_1,ind_node_3) = 1;
A(ind_node_3,ind_node_1) = 1;

A(ind_node_2,ind_node_4) = 1;
A(ind_node_4,ind_node_2) = 1;
end

function B = compute_B(A,C,N_comm,IN) % Compute a modified incidence matrix B
N = size(A,1);
Deg = A*ones(N,1);
B = [];
for counter_1 = 1 : N - 1
    for counter_2 = counter_1 + 1 : N
        if(A(counter_1,counter_2))
            B = [B, [counter_1;counter_2;IN(counter_1);IN(counter_2);...
                 C(counter_1);C(counter_2);Deg(counter_1);Deg(counter_2);...
                 N_comm(C(counter_1));N_comm(C(counter_2))]];
        end
    end
end
end

function Q = compute_modularity(A,C) % Compute modularity of the given partition
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*C - C'*ones(1,N)) == 0)));
end
