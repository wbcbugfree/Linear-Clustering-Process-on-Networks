function [C,c,m] = LCP_c(A,c)
% This function implements the Linear Clustering Process (LCP), from the
% publication
N = size(A,1);                                                                          % Number of nodes
N_pr = 40;                                                                              % Number of iterations
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
    A_model = eye(N) + A_fm.*W - diag(A_fm.*W*ones(N,1));                               % Compute state space matrix                                           % Compute the position vector in the steady state
    [ind_r,pos_r,~] = compute_y_2_v2(A_model);                                          % Compute the position vector in the steady state
    Inds(:,counter) = ind_r; Pos_M(:,counter) = pos_r;                                  % Store the position of each node per iteration
    [~,m_it(counter),C_it(:,counter),c_it(counter)] = ident_communities_v1_N(A(ind_r,ind_r),c);   % Identify the clusters, for given position vector
    [A_fm,~,~] = scale_links(A_fm,ind_r,L_rem(counter),counter,N_pr);                   % Scale the weights of the identified inter-community links links and recompute the adjacency matrix
    counter = counter + 1;                                                              % Update the counter
end
[val_mod,ind] = max(m_it);                                                              % Adopt the iteration with the highest modularity
c = c_it(ind);                                                                          % Output the estimated number of clusters
m = val_mod;                                                                            % Output the estimated modularity
C = zeros(N,1);
C(Inds(:,ind)) = C_it(:,ind);                                                           % Output the cluster membership function of each node                                                                               % Return C and m
end

%% Used sub-functions
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
W = (alpha + delta)*Deg_inv*(A_s + A)*Deg_inv - 0.5*delta*(Deg_inv*A + A*Deg_inv);      % Compute matrix W (Theorem 1, page 6)
% W = (alpha + 1/(1+alpha)^4)*Deg_inv*(A_s + A)*Deg_inv - 0.5*1/(1+alpha)^4*(Deg_inv*A + A*Deg_inv);      % Compute matrix W (Theorem 1, page 6)
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

function [A_new,A_rem,Deg_new] = scale_links(A,ind,trsh,br_it,N_pr)
% Adjacency matrix A provided as input is not relabeled!
N = size(A,1);                                                      % Network size
A_n = A(ind,ind);                                                   % Relabeled adjacency matrix A based on y_2 (eq. 10, page 11)
Deg_v = A_n*ones(N,1);                                              % Relabeled degree vector based on y_2 (eq. 10, page 11)
N_rem_1 = ceil(nnz(A)/2*trsh/100);                                  % Number of links, whose weight to scale wihtin one iteration
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

function [M,m,C_out,N_clusters] = ident_communities_v1_N(A,N_cl)
% Estimate communities in the network
N = size(A,1);
Deg = A*ones(N,1);
% Step 1: Find cluster borders
L_2 = nnz(A);
C = ident_clusters_v2_N(A,N,Deg,L_2,1,ceil(log2(N_cl))+1);
if(isempty(C))
    M = 0;m = 0;C_out=ones(N,1);N_clusters = 1;return
end
A_mat = A - (1/nnz(A)).*(Deg*Deg');
%% Step 2: Compute the modularity matrix M
[M,C_out,N_clusters] = compute_modularity_matrix_v1(C,A_mat,A,N);
%% Step 3: Optimise partition by maximising modularity
[M,C_out,N_clusters] = optimise_clusters_v1_N(M,C_out,N_clusters,N_cl);
% m = trace(M)
m = compute_modularity_v1(A,C_out);
end

function C = ident_clusters_v2_N(A,N,Deg,L_2,br_it,N_it)
Mod_1 = zeros(N,1);         Mod_2 = zeros(N,1);
Mod_1(1) = -Deg(1)^2/(L_2);   Mod_2(N) = -Deg(N)^2/(L_2);
A_mat = A - (Deg*Deg')./L_2;
for br = 2 : N
    Mod_1(br) = Mod_1(br-1) + 2*sum(A_mat(1:br-1,br)) + A_mat(br,br);
    Mod_2(N-br+1) = Mod_2(N-br+2) + 2*sum(A_mat(N-br+2:N,N-br+1)) + A_mat(N-br+1,N-br+1);
end
[m_max,ind_max] = max(Mod_1+Mod_2);       % Determine cluster border
m_max = m_max(1); ind_max = ind_max(1);
if((ind_max == 1 || ind_max == N) && br_it < N_it && N > 2)
    ind_max = ceil(N/2);
    A_1 = A(1:ind_max,1:ind_max); Deg_1 = Deg(1:ind_max); N_1 = ind_max;
    A_2 = A(ind_max+1:N,ind_max+1:N); Deg_2 = Deg(ind_max+1:N); N_2 = N - ind_max;
    C = [ident_clusters_v2_N(A_1,N_1,Deg_1,L_2,br_it + 1,N_it),ind_max,...
         ind_max + ident_clusters_v2_N(A_2,N_2,Deg_2,L_2,br_it + 1,N_it)];
elseif(ind_max > 1 && ind_max < N && br_it < N_it && N > 2)% Is there a new cluster border?
    A_1 = A(1:ind_max,1:ind_max); Deg_1 = Deg(1:ind_max); N_1 = ind_max;
    A_2 = A(ind_max+1:N,ind_max+1:N); Deg_2 = Deg(ind_max+1:N); N_2 = N - ind_max;
    C = [ident_clusters_v2_N(A_1,N_1,Deg_1,L_2,br_it + 1,N_it),ind_max,...
         ind_max + ident_clusters_v2_N(A_2,N_2,Deg_2,L_2,br_it + 1,N_it)];
else
    C = [];
end;end

function [M,C_vec,N_cl] = compute_modularity_matrix_v1(C,A_mat,A,N)
Inv = 1/nnz(A);     % 1/2L
N_cl = length(C)+1; % Number of clusters
C_low = [1,C]; C_upp = [C,N]; C_vec = zeros(N,1); 
for br = 1 : N_cl   % Construct cluster membership vector C
    C_vec(C_low(br):C_upp(br)) = br;
end
M = zeros(N_cl);    % Construct the cluster membership matrix M
for br_1 = 1 : N_cl
    for br_2 = 1 : br_1
        M(br_1,br_2) = Inv*sum(sum(A_mat(C_low(br_1):C_upp(br_1),C_low(br_2):C_upp(br_2))));
        M(br_2,br_1) = M(br_1,br_2);
    end
end
end

function [M,C_vec,N_cl] = optimise_clusters_v1_N(M,C_vec,N_cl,c)
Test = diag(M,1);
while(length(Test) > c-1)                                           % Find position of the largest positive off diagonal element of M
    ind = find(Test == max(Test));
    for br = ind + 1 : N_cl                                         % Update the cluster membership vector
        C_vec(C_vec == br) = br - 1;
    end
    N_cl = N_cl - 1;                                                % Update the number of clusters
    M(:,ind) = M(:,ind) + M(:,ind+1); M(:,ind+1) = [];              % Update the cluster membership matrix M
    M(ind,:) = M(ind,:) + M(ind+1,:); M(ind+1,:) = [];
    Test = diag(M,1);                                               % Update the second diagonal vector
end;end

function Q = compute_modularity_v1(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end
