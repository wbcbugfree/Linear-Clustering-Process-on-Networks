function [C,c,m] = LCP_v2(A)
% This function implements the Linear Clustering Process (LCP), from the
% publication
N = size(A,1);                                                                          % Number of nodes
N_pr = 30;                                                                              % Number of iterations
L_per = 60/N_pr;                                                                        % Ratio of links whose weight will be scaled
L_rem = linspace(L_per,L_per,N_pr);                                                     % Scaling intensity for links in each iteration

alpha = 0.95;                                                                           % Attraction strength
delta = 1e-3;                                                                           % Repulsion strength

A_fm = A;                                                                               % Store the adjacency matrix

Inds = zeros(N,N_pr);                                                                   % Initalise the ranking of each node in y_2 per iteration
Pos_M = zeros(N,N_pr);                                                                  % Initalise the position of each node in y_2 per iteration
m_it = zeros(N_pr,1);                                                                   % Initalise the modularity in each iteration
C_it = zeros(N,N_pr);                                                                   % Initalise the number of clusters in each iteration
c_it = zeros(N_pr,1);                                                                   % Initalise the modularity in each iteration

W = compute_W(A_fm,N,alpha,delta);                                                      % Compute W matrix (we compute it only once!)
counter = 1;
while(counter <= N_pr)                                                                  % Perform N_pr iterations of LCP
    A_model = eye(N) + A_fm.*W - diag(A_fm.*W*ones(N,1));                               % Compute state space matrix
    [ind_r,pos_r,y_2] = compute_y_2(A_model);                                           % Compute the position vector in the steady state
    %     [ind_r,pos_r,y_2] = compute_y_2_v2(A_model);                                  % Compute the position vector in the steady state
    Inds(:,counter) = ind_r; Pos_M(:,counter) = pos_r;                                  % Store the position of each node per iteration
    [~,m_it(counter),C_it(:,counter),c_it(counter)] = optimise_m(A(ind_r,ind_r));       % Identify the clusters, for given position vector
    %     [A_fm,~,~] = scale_links(A_fm,ind_r,L_rem(counter),counter,N_pr,y_2);         % Scale the weights of the identified inter-community links links and recompute the adjacency matrix
    C = ones(N,1)*C_it(:,counter)' - C_it(:,counter)*ones(1,N); C = (abs(C) > 0);
    [A_fm,~,~] = scale_links_v2(A_fm,ind_r,L_rem(counter),counter,N_pr,y_2,C);          % Scale the weights of the identified inter-community links links and recompute the adjacency matrix
    counter = counter + 1;                                                              % Update the counter
end
[val_mod,ind] = max(m_it);                                                              % Adopt the iteration with the highest modularity
c = c_it(ind);                                                                          % Output the estimated number of clusters
m = val_mod;                                                                            % Output the estimated modularity
C = zeros(N,1);
C(Inds(:,ind)) = C_it(:,ind);                                                           % Output the cluster membership function of each node

% Plot the results per iteration
%{
figure
subplot(1,2,1)
plot(m_it)
xlabel('Iteration','interpreter','latex')
ylabel('$m$','interpreter','latex')
title('Estimated modularity $m$ per iteration','interpreter','latex')
set(gca,'Fontsize',14,'TickLabelInterpreter','latex')
subplot(1,2,2)
plot(c_it)
xlabel('Iteration','interpreter','latex')
ylabel('$c$','interpreter','latex')
title('Estimated number of clusters $c$ per iteration','interpreter','latex')
set(gca,'Fontsize',14,'TickLabelInterpreter','latex')
%}
return                                                                                  % Return C and m
end

%% Used sub-functions
function [A_new,A_rem,Deg_new] = scale_links(A,ind,trsh,br_it,N_pr,y_2)
% Adjacency matrix A provided as input is not relabeled!
N = size(A,1);                                                      % Network size
A_n = A(ind,ind);                                                   % Relabeled adjacency matrix A based on y_2 (eq. 10, page 11)
Deg_v = A_n*ones(N,1);                                              % Relabeled degree vector based on y_2 (eq. 10, page 11)
N_rem_1 = ceil(nnz(A)/2*trsh/100);                                  % Number of links, whose weight to scale wihtin one iteration
% Dist_M_1 = (abs(diag(1:N)*A_n - A_n*diag(1:N))).*A_n;               % Compute the distance matrix between any two adjacent nodes (sec 5.1, page 15)
Dist_M_1 = (abs(diag(1:N)*A_n - A_n*diag(1:N))).*A_n;               % Compute the distance matrix between any two adjacent nodes (sec 5.1, page 15)
A_new = A_n;                                                        % Store the adjacency matrix
counter_1 = 1;                                                      % Coutner for the number of links with scaled weight
while counter_1 < N_rem_1                                           % Run the loop until N_rem links are scaled in weight
    tr = 1;                                                         % Initialise the indicator function in case the link is found.
    while (tr)                                                      % Run the loop until the link is removed
        val = max(Dist_M_1(:));                                     % Find the maximum distance value
        [ind_1,ind_2] = find(Dist_M_1 == val);                      % Identify the link connecting two most far away nodes
        if((Deg_v(ind_1(1))) > 1 && (Deg_v(ind_2(1))) > 1 ...
                && A_new(ind_1(1),ind_2(1)) == 1)                       % Check the validity of a link (i.e. if adjacent nodes have degree > 1)?
            Dist_M_1(ind_1(1),ind_2(1)) = 0;                        % Define the distance of a chosen link as 0, so it is not conisdered anymore
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
            counter_1 = counter_1 + 1; tr = 0;                      % Increment the counter
            Deg_v(ind_1(1)) = Deg_v(ind_1(1)) - 1;                  % Update the degree values for each node
            Deg_v(ind_2(1)) = Deg_v(ind_2(1)) - 1;
            A_new(ind_1(1),ind_2(1)) = (0.01/(N_pr))*br_it;        % Scale the weight of a link
            A_new(ind_2(1),ind_1(1)) = (0.01/(N_pr))*br_it;
            A_new(ind_1(1),ind_2(1)) = 0;        % Scale the weight of a link
            A_new(ind_2(1),ind_1(1)) = 0;
        else                                                        % This link can not be removed
            Dist_M_1(ind_1(1),ind_2(1)) = 0;                        % Define the distance of a chosen link as 0, so it is not conisdered anymore
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
        end
    end
end
A_new(ind,ind) = A_new;
A_rem = (A - A_new);
Deg_new = A_new*ones(N,1);
% Check if there are zero-degree nodes (testing the code)
if(nnz(Deg_new) < N)
    disp('Number of zero degree nodes:')
    disp(length(Deg_new) - nnz(Deg_new))
    plot(sort(Deg_new))
    error('Something is wrong!!!')
end
end

function [A_new,A_rem,Deg_new] = scale_links_v2(A,ind,trsh,br_it,N_pr,y_2,C)
% Adjacency matrix A provided as input is not relabeled!
N = size(A,1);                                                      % Network size
A_n = A(ind,ind);                                                   % Relabeled adjacency matrix A based on y_2 (eq. 10, page 11)
Deg_v = A_n*ones(N,1);                                              % Relabeled degree vector based on y_2 (eq. 10, page 11)
N_rem_1 = ceil(nnz(A)/2*trsh/100);                                  % Number of links, whose weight to scale wihtin one iteration
Dist_M_1 = (abs(diag(1:N)*A_n - A_n*diag(1:N))).*A_n;             % Compute the distance matrix between any two adjacent nodes (sec 5.1, page 15)
% Dist_M_1 = (ones(N) + abs(diag(1:N)*A_n - A_n*diag(1:N))).*A_n.*C;  % Compute the distance matrix between any two adjacent nodes (sec 5.1, page 15)
A_new = A_n;                                                        % Store the adjacency matrix
counter_1 = 1;                                                      % Coutner for the number of links with scaled weight
while counter_1 < N_rem_1                                           % Run the loop until N_rem links are scaled in weight
    tr = 1;                                                         % Initialise the indicator function in case the link is found.
    while (tr)                                                      % Run the loop until the link is removed
        Dist_v = Dist_M_1(:); Dist_v = Dist_v(Dist_v > 0);
        val = max(Dist_v);                                     % Find the maximum distance value
        %         val = max(Dist_M_1(:));                                     % Find the maximum distance value
        [ind_1,ind_2] = find(Dist_M_1 == val(1));                      % Identify the link connecting two most far away nodes
        if((Deg_v(ind_1(1))) > 1 && (Deg_v(ind_2(1))) > 1 ...
                && A_new(ind_1(1),ind_2(1)) == 1)                       % Check the validity of a link (i.e. if adjacent nodes have degree > 1)?
            Dist_M_1(ind_1(1),ind_2(1)) = 0;                        % Define the distance of a chosen link as 0, so it is not conisdered anymore
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
            counter_1 = counter_1 + 1; tr = 0;                      % Increment the counter
            Deg_v(ind_1(1)) = Deg_v(ind_1(1)) - 1;                  % Update the degree values for each node
            Deg_v(ind_2(1)) = Deg_v(ind_2(1)) - 1;
            A_new(ind_1(1),ind_2(1)) = (0.001/(N_pr))*br_it;        % Scale the weight of a link
            A_new(ind_2(1),ind_1(1)) = (0.001/(N_pr))*br_it;
        else                                                        % This link can not be removed
            Dist_M_1(ind_1(1),ind_2(1)) = 0;                        % Define the distance of a chosen link as 0, so it is not conisdered anymore
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
        end
    end
end
A_new(ind,ind) = A_new;
A_rem = (A - A_new);
Deg_new = A_new*ones(N,1);
% Check if there are zero-degree nodes (testing the code)
if(nnz(Deg_new) < N)
    disp('Number of zero degree nodes:')
    disp(length(Deg_new) - nnz(Deg_new))
    plot(sort(Deg_new))
    error('Something is wrong!!!')
end
end

function [ind_rank,pos_rank,y_2] = compute_y_2(W)
[Y,B] = eig(W);                     % Compute the eigenvalue decomposition of I + W - diag(W*u)
[~,ind_pe] = maxk(diag(B),3);       % Identify the three largest eigenvalues
y_2 = Y(:,ind_pe(2));               % Store the eigenvector y_2
[pos_rank,ind_rank] = sort(y_2);    % Compute position and ranking of eacn node in y_2
end

function W = compute_W(A,N,alpha,delta)
A_s = zeros(N);                                                                         % Initialise A.*A^2
for i = 1 : N                                                                           % Compute A.*A^2 as in algorithm 2 (page 17)
    ind_n = find(A(i,:) == 1);                                                          % (pseudocode 2, line 3)
    for j = 1  : length(ind_n)                                                          % (pseudocode 2, line 3)
        ind_n_2 = find(A(ind_n(j),:) == 1); ind_p_2 = ind_n_2(ind_n_2 > i);             % (pseudocode 2, line 4)
        A_s(i,ind_p_2) = A_s(i,ind_p_2) + A(i,ind_p_2).*A(ind_n(j),ind_p_2);            % (pseudocode 2, line 5)
    end
end
A_s = A_s + A_s';                                                                       % Compute A.*A^2 (pseudocode 2, line 9)
Deg_inv = diag(A*ones(N,1))^-1;                                                         % Compute vector with inverse degrees of each node
W = (alpha + delta)*Deg_inv*(A_s + A)*Deg_inv + 0.5*delta*(Deg_inv*A + A*Deg_inv);      % Compute matrix W (Theorem 1, page 6)
end

function [M,m,C_out,N_clusters] = optimise_m(A)
% This function computes merger matrix M, number of clusters c and cluster
% membership function C
N = size(A,1);                                                      % Number of nodes N
Deg = A*ones(N,1);                                                  % Degree vector
L_2 = nnz(A);                                                       % Twice number of links
C = ident_clusters_border(A,N,Deg,L_2,0);                           % Identify cluster membership of each node
if(isempty(C))                                                      % If there are no clusters
    M = 0;m = 0;C_out=ones(N,1);N_clusters = 1;return               % Return only one cluster, with all nodes in it
end
A_mat = A - (1/nnz(A)).*(Deg*Deg');                                 % Compute the modularity matrix A - 1/2L*d*d^T
[M,C_out,N_clusters] = compute_M(C,A_mat,A,N);                      % Compute the modularity matrix M
[M,C_out,N_clusters] = merge_clusters(M,C_out,N_clusters);          % Optimise partition by maximising modularity
m = compute_modularity_m(A,C_out);                                  % Compute the modularity for a given partition
end

function C = ident_clusters_border(A,N,Deg,L_2,trsh)
% This function implements the algorithm 1 from the LCP paper (page 13)
Mod_1 = zeros(N,1);         Mod_2 = zeros(N,1);                                             % Initialise the modularity vectors (pseudocode 1, line 2)
Mod_1(1) = -Deg(1)^2/(L_2); Mod_2(N) = -Deg(N)^2/(L_2);                                     % Store the first value of the modularity vectors (pseudocode 1, line 4-5)
A_mat = A - (Deg*Deg')./L_2;
for br = 2 : N                                                                              % Iteratively compute the modularity of each possible partition (pseudocode 1, lines 7 - 15)
    Mod_1(br) = Mod_1(br-1) + 2*sum(A_mat(1:br-1,br)) + A_mat(br,br);                       % Update the modularity of each possible bisection (pseudocode 1, line 9,11,13)
    Mod_2(N-br+1) = Mod_2(N-br+2) + 2*sum(A_mat(N-br+2:N,N-br+1)) + A_mat(N-br+1,N-br+1);   % Update the modularity of each possible bisection (pseudocode 1, line 10,12,14)
end
[m_max,ind_max] = max(Mod_1+Mod_2);                                                         % Determine the best partition into two clusters       (pseudocode 1, line 16)
m_max = m_max(1); ind_max = ind_max(1);                                                     % Store the cluster border and the obtained modularity (pseudocode 1, line 16)
if(m_max >= trsh && ind_max > 1 && ind_max < N)                                             % If the new partition improves modularity of the parent cluster, adopt it (pseudocode 1, line 17)
    A_1 = A(1:ind_max,1:ind_max); Deg_1 = Deg(1:ind_max); N_1 = ind_max;                    % Perform the partition and compute adjacency matrix and degree vector for the first cluster (pseudocode 1, line 18)
    A_2 = A(ind_max+1:N,ind_max+1:N); Deg_2 = Deg(ind_max+1:N); N_2 = N - ind_max;          % Perform the partition and compute adjacency matrix and degree vector for the second cluster (pseudocode 1, line 18)
    C = [ident_clusters_border(A_1,N_1,Deg_1,L_2,Mod_1(ind_max)),ind_max,...
        ind_max + ident_clusters_border(A_2,N_2,Deg_2,L_2,Mod_2(ind_max))];                % Call the function to check if the obtained two clusters can be further partitioned (pseudocode 1, line 20)
else                                                                                        % otherwise, adopt only the parent cluster (pseudocode 1, line 21)
    C = [];
end
end

function [M,C_vec,N_cl] = compute_M(C,A_mat,A,N)
% This function computes the c x c modularity matrix M and the Nx1 cluster
% membership vector C
Inv = 1/nnz(A);                                     % 1/2L
N_cl = length(C)+1;                                 % Number of clusters
C_low = [1,C]; C_upp = [C,N]; C_vec = zeros(N,1);   % Explanation
for br = 1 : N_cl                                   % Construct cluster membership vector C
    C_vec(C_low(br):C_upp(br)) = br;
end
M = zeros(N_cl);                                    % Construct the cluster membership matrix M
for br_1 = 1 : N_cl
    for br_2 = 1 : br_1
        M(br_1,br_2) = Inv*sum(sum(A_mat(C_low(br_1):C_upp(br_1),C_low(br_2):C_upp(br_2))));
        M(br_2,br_1) = M(br_1,br_2);
    end
end
end

function [M,C_vec,N_cl] = merge_clusters(M,C_vec,N_cl)
% This function further optimises the modularity m
Off_diag_M = diag(M,1);                                             % Store the elements above the main diagonal of the modularity matrix M
while((max(Off_diag_M) > 0))                                        % While there is a positive off-diagonal element in M
    ind = find(Off_diag_M == max(Off_diag_M));                      % Find position of the largest off-diagonal element in M
    for br = ind + 1 : N_cl                                         % Since we merge two clusters, update memebership of all nodes affected
        C_vec(C_vec == br) = br - 1;                                % Update the cluster membership vector C
    end
    N_cl = N_cl - 1;                                                % Update number of clusters
    M(:,ind) = M(:,ind) + M(:,ind+1); M(:,ind+1) = [];              % Update the cluster membership matrix M
    M(ind,:) = M(ind,:) + M(ind+1,:); M(ind+1,:) = [];              % Update the cluster membership matrix M
    Off_diag_M = diag(M,1);                                         % Update the off-diagonal elements of M
end
end

function Q = compute_modularity_m(A,C)
N = size(A,1); Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*C' - C*ones(1,N)) == 0)));
end

function [ind_rank,pos_rank,y_2] = compute_y_2_v2(A_model)
N = size(A_model,1);
A_m = A_model - 1/N*ones(N);
Vec = rand(N,1);
ind = 1;
cnt = 1;
while(ind)
    cnt = cnt + 1;
    Vec_est = A_m*Vec;
    Vec_est = Vec_est./(sqrt(sum(Vec_est.*Vec_est)));
    if(mean(abs(Vec - Vec_est)) < 1e-5)
        ind = 0;
    else
        Vec = Vec_est;
    end
end
[pos_rank,ind_rank] = sort(Vec);    % Compute position and ranking of eacn node in y_2
y_2 = Vec;
end