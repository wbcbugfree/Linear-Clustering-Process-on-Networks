% This code runs the proposed iterative linear process and estimates the
% partition
clc
clear
close all

%% Step 1: Define the SBM network parameters
N_clusters = 8;
N = 1000;                                                               % Number of nodes
d_av = 7;

% b_in = linspace(8,14,20);   % c = 2 clusters
% b_in = linspace(9,20,20);   % c = 3 clusters
% b_in = linspace(11,25,20);   % c = 4 clusters
b_in = linspace(21,55,20);   % c = 8 clusters
% b_in = linspace(17,62,20);   % c = 10 clusters
% b_in = linspace(45,125,20);   % c = 20 clusters

b_out = (d_av*N_clusters - b_in)./(N_clusters - 1);

% b_in - b_out
% b_out
% pause()

N_it = 1e3;

N_LCP = zeros(length(b_in),N_it);
M_LCP = zeros(length(b_in),N_it);

N_LCP_N = zeros(length(b_in),N_it);
M_LCP_N = zeros(length(b_in),N_it);

N_LCP_B = zeros(length(b_in),N_it);

N_louv = zeros(length(b_in),N_it);
M_louv = zeros(length(b_in),N_it);

N_newm = zeros(length(b_in),N_it);
M_newm = zeros(length(b_in),N_it);

N_nbt = zeros(length(b_in),N_it);

M_orig = zeros(length(b_in),N_it);


for br = 1 : length(b_in)
    for br_2 = 1 : N_it
        [G,G_bidir,links,links_dir,nodes,groups] = generate_sbm(1000,N_clusters,b_in(br),b_out(br));
        A = adjacency(G);
        N = size(A,1);
        Deg = A*ones(N,1);

        ind = find(Deg>0);
        A = A(ind,ind);
        N = length(ind);
        Deg = Deg(ind);

        C = groups(ind)';
        [val,pos] = sort(C);
        C = C(pos);
        A_0 = A(pos,pos);
        Q_0 = compute_modularity(A_0,C);
        M_orig(br,br_2) = Q_0;

        [N_LCP_N(br,br_2),M_LCP_N(br,br_2),N_LCP(br,br_2),M_LCP(br,br_2)] = linear_clustering_orig_and_N(A,N_clusters);
        
        N_LCP_B(br,br_2) = Non_back_tracking_LCP(A,N,0.95);

        Res = cluster_jl(A);
        val_mod = Res.MOD; [val_louv,ind_louv] = max(val_mod);

        M_louv(br,br_2) = val_louv;
        N_louv(br,br_2) = length(Res.SIZE{ind_louv(1)});
%         M_louv(br,br_2) = Res.MOD(end);
%         N_louv(br,br_2) = length(Res.SIZE{end});

        [N_newm(br,br_2),M_newm(br,br_2)] = Newman_clustering(A);
        N_nbt(br,br_2) = Non_back_tracking(links_dir,A,N);

    end
end

N_LCP_mean = zeros(20,1);
N_LCP_B_mean = zeros(20,1);
N_louv_mean = zeros(20,1);
N_newm_mean = zeros(20,1);
N_nbt_mean = zeros(20,1);

M_orig_mean = zeros(20,1);
M_LCP_mean = zeros(20,1);
M_LCP_N_mean = zeros(20,1);
M_louv_mean = zeros(20,1);
M_newm_mean = zeros(20,1);

for br = 1 : 20
    N_LCP_mean(br)      = mean(N_LCP(br,:));
    N_LCP_B_mean(br)    = mean(N_LCP_B(br,:));
    N_louv_mean(br)     = mean(N_louv(br,:));
    N_newm_mean(br)     = mean(N_newm(br,:));
    N_nbt_mean(br)      = mean(N_nbt(br,:));

    M_orig_mean(br)     = mean(M_orig(br,:));
    M_LCP_mean(br)      = mean(M_LCP(br,:));
    M_LCP_N_mean(br)    = mean(M_LCP_N(br,:));
    M_louv_mean(br)     = mean(M_louv(br,:));
    M_newm_mean(br)     = mean(M_newm(br,:));
end

N_est_mean = [N_LCP_mean, N_LCP_B_mean, N_louv_mean, N_newm_mean, N_nbt_mean];
M_est_mean = [M_orig_mean, M_LCP_mean, M_LCP_N_mean, M_louv_mean, M_newm_mean];

% N_est_All = [mean(N_LCP') mean(N_LCP_B') mean(N_louv') mean(N_newm') mean(N_nbt')];
% M_est_All = [mean(M_LCP') mean(M_LCP_N') mean(M_louv') mean(M_newm')];

disp('Estimated Number of Clusters')
disp(N_est_mean)
disp('Estimated modularity')
disp(M_est_mean)
% figure
% plot(b_in-b_out,N_clusters*ones(length(b_in),1),'k','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(N_LCP'),'b-x','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(N_louv'),'r-+','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(N_newm'),'g-*','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(N_LCP_B'),'m-<','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(N_nbt'),'c-d','LineWidth',1.1)
% xline(sqrt(d_av)*N_clusters,'--k',{'Detectability','Limit'});
% xlim([0.91*min(b_in-b_out),1.02*max(b_in-b_out)])
% legend('$c$','LCP','Louvain','Newman','LCP$_{h}$','Non-back tracking','interpreter','latex','Color','none','box','off','Location','northeast')
% xlabel('$b_{in}-b_{out}$','interpreter','latex')
% ylabel('Number of clusters','interpreter','latex')
% title('Estimated number of clusters','interpreter','latex')
% set(gca,'Fontsize',14,'TickLabelInterpreter','latex')

% figure
% plot(b_in-b_out,mean(M_orig'),'k','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(M_LCP'),'b-x','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(M_louv'),'r-+','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(M_newm'),'g-*','LineWidth',1.1)
% hold on
% plot(b_in-b_out,mean(M_LCP_N'),'m-<','LineWidth',1.1)
% xline(sqrt(d_av)*N_clusters,'--k',{'Detectability','Limit'});
% xlim([0.91*min(b_in-b_out),1.02*max(b_in-b_out)])
% legend('$m$','LCP','Louvain','Newman','LCP$_{n}$','interpreter','latex','Color','none','box','off','Location','northwest')
% xlabel('$b_{in}-b_{out}$','interpreter','latex')
% ylabel('Modularity','interpreter','latex')
% title('Modularity of the estimated partition','interpreter','latex')
% set(gca,'Fontsize',14,'TickLabelInterpreter','latex')

%% Used functions

function [G,G_bidir,links,links_dir,nodes,groups] = generate_sbm(n,c,b_in,b_out)
%GENERATE_SBM Generate SBM graph
%   Generate a graph with n nodes according to the Symmetric Stochastic
%   Block Model with c groups

b = b_out + ((b_in - b_out)/c);  % expected average degree
sizes = fix(n/c) + zeros(c,1);  % equal parts groups
sizes(1:mod(n,c)) = sizes(1:mod(n,c)) + 1;  % accomodate for imbalanced groups
nodes = 1:n;
groups_dist = repelem(1:c,sizes);  % create vector with groups according to sizes
groups = groups_dist(randperm(n));  % randomly distribute groups

% Generate links
links = [];
for i = 1:(n-1)
    for j = (i+1):n  % starting from (i+1) because no self-loops
        has_link = 0;
        if groups(i) == groups(j)
            if rand() < (b_in/n)
                has_link = 1;
            end
        else
            if rand() < (b_out/n)
                has_link = 1;
            end
        end
        if has_link
            links = [links; [i j]];
        end
    end
end
links_dir = sortrows([links(:,1) links(:,2); links(:,2) links(:,1)]);

% Random labeling

G = graph(links(:,1),links(:,2));
G_bidir = digraph(links_dir(:,1),links_dir(:,2));

end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end

function [C_N,m_N,C_orig,m_orig] = linear_clustering_orig_and_N(A,c)

% This function iteratively implements the theoretical solution of the proposed
% model, in order to recover clusters in a network
N = size(A,1);
N_pr = 31;
L_per = 60/N_pr;
L_rem = linspace(L_per,L_per,N_pr);

alpha = 0.98;
delta = 1e-4;

A_fm = A;
Deg_fm = A_fm*ones(N,1);

Inds = zeros(N,N_pr);
Pos_M = zeros(N,N_pr);

m_new_N = zeros(N_pr,1);
C_new_N = zeros(N,N_pr);
N_c_new_N = zeros(N_pr,1);

m_new_orig = zeros(N_pr,1);
C_new_orig = zeros(N,N_pr);
N_c_new_orig = zeros(N_pr,1);

% Compute W
W = compute_A_model_v4(A_fm,Deg_fm,N,alpha,delta);

br = 1;
while(br <= N_pr)
    % Compute state space matrix
    A_model = eye(N) + A_fm.*W - diag(A_fm.*W*ones(N,1));
    % Compute the position vector in the steady state
    [ind_r,pos_r,~] = compute_ranking_f_1(A_model,A);
    % Store the rankings of nodes
    Inds(:,br) = ind_r;
    Pos_M(:,br) = pos_r;
    % Identify the clusters, for given position vector
    [~,m_new_orig(br),C_new_orig(:,br),N_c_new_orig(br)] = ident_communities_v3(A(ind_r,ind_r));
    [~,m_new_N(br),C_new_N(:,br),N_c_new_N(br)] = ident_communities_v1_N(A(ind_r,ind_r),c);

    % Remove the links and recompute the solution
    [A_new,~,~] = remove_links_6(A_fm,ind_r,L_rem(br),br,N_pr);
    % Update the adjacency matrix
    A_fm = A_new;  
    br = br + 1;
end

[val_mod_N,ind_N] = max(m_new_N);
C_N = N_c_new_N(ind_N);
m_N = val_mod_N;

[val_mod_orig,ind_orig] = max(m_new_orig);
C_orig = N_c_new_orig(ind_orig);
m_orig = val_mod_orig;
end

function [A_new,A_rem,Deg_new] = remove_links_6(A,ind,trsh,br_it,N_pr,M,C)
% Adjacency matrix A provided as input is not relabeled.
N = size(A,1);
A_n = A(ind,ind);

Deg_v = A_n*ones(N,1);
N_rem_1 = ceil(nnz(A)/2*trsh/100);

Dist_M_1 = (abs(diag(1:N)*A_n - A_n*diag(1:N))).*A_n;
A_new = A_n;

brojac_1 = 1;
while brojac_1 < N_rem_1
    tr = 0;
    br_2 = 0;
    % Run the loop untile the link is removed
    while (tr == 0)
        br_2 = br_2 + 1;
        val = max(Dist_M_1(:));
        [ind_1,ind_2] = find(Dist_M_1 == val);
        if((Deg_v(ind_1(1))) > 1 && (Deg_v(ind_2(1))) > 1 && A_new(ind_1(1),ind_2(1)) == 1)
            % Remove the link
            Dist_M_1(ind_1(1),ind_2(1)) = 0;
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
            brojac_1 = brojac_1 + 1; tr = 1;
            Deg_v(ind_1(1)) = Deg_v(ind_1(1)) - 1;
            Deg_v(ind_2(1)) = Deg_v(ind_2(1)) - 1;
            % Scale the weights
            A_new(ind_1(1),ind_2(1)) = (0.05/(N_pr))*br_it;
            A_new(ind_2(1),ind_1(1)) = (0.05/(N_pr))*br_it;
        else
            % This link can not be removed
            Dist_M_1(ind_1(1),ind_2(1)) = 0;
            Dist_M_1(ind_2(1),ind_1(1)) = 0;
        end
    end
end

A_new(ind,ind) = A_new;
A_rem = (A - A_new);
Deg_new = A_new*ones(N,1);

if(nnz(Deg_new) < N)
    disp('Number of zero degree nodes:')
    disp(length(Deg_new) - nnz(Deg_new))
    plot(sort(Deg_new))
    error('Something is wrong!!!')
end
end

function [ind_rank,pos_rank,pos_v] = compute_ranking_f_1(A_model,A)

[X_mod,L_mod] = eig(A_model);
% Eig_v = abs(diag(L_mod) - 1);
% [~,ind_pe] = mink(Eig_v,20);
[~,ind_pe] = maxk(diag(L_mod),3);

% figure
% % plot(sort(X_mod(:,ind_pe(1))))
% pause()
% close
pos_v = X_mod(:,ind_pe(2));% + X_mod(:,ind_pe(3));
% pos_v = X_mod(:,ind_pe(2)) + (L_mod(ind_pe(3),ind_pe(3))/L_mod(ind_pe(2),ind_pe(2)))*X_mod(:,ind_pe(3));

% pos_v = 1/L_mod(ind_pe(2),ind_pe(2))*X_mod(:,ind_pe(2)) + 1/L_mod(ind_pe(3),ind_pe(3))*X_mod(:,ind_pe(3));
% [pos_rank,ind_rank] = sort( 1*X_mod(:,ind_pe(1)) + 0*X_mod(:,ind_pe(2)));
[pos_rank,ind_rank] = sort(pos_v);
% Q_rk = compute_modularity(A(ind_rank,ind_rank),C);

% N = size(A,1);
% [X_mod,~] = eigs(A_model-1/N*ones(N),10);
% pos_v = X_mod(:,2);
% [pos_rank,ind_rank] = sort(pos_v);
% Q_rk = compute_modularity(A(ind_rank,ind_rank),C);

% N = size(A,1);
% Vec_st = rand(N,1);
% [vec,~]=power_method(Vec_st/(Vec_st'*Vec_st),A_model - (1/N)*ones(N),1e-8);
% % [ vec, ~, ~ ] = power_method_3 ( N, A_model - 1/N*ones(N), rand(N,1), 1000, 1e-12);
% % [ vec, ~, ~ ] = power_method_3 ( N, A_model, N*rand(N,1), 800, 1e-10);
% [pos_rank,ind_rank] = sort(vec);
% pos_v = vec;
% Q_rk = compute_modularity(A(ind_rank,ind_rank),C);
end

function W = compute_A_model_v4(A,~,N,alpha,delta)
% A_s = (A^2);
A_s = zeros(N);
for i = 1 : N
    ind_n = find(A(i,:) == 1);
    for j = 1  : length(ind_n)
        ind_n_2 = find(A(ind_n(j),:) == 1);
        ind_p_2 = ind_n_2(ind_n_2 > i);
        A_s(i,ind_p_2) = A_s(i,ind_p_2) + A(i,ind_p_2).*A(ind_n(j),ind_p_2);
    end
end
% A_s = A_s + A_s';
A_s = zeros(N);
for br_i = 1 : N
    ind_n = find(A(br_i,1:N) == 1);
    A_s(ind_n,ind_n) = A_s(ind_n,ind_n) + 1;
end
V_s = diag(A_s);

T_a_att = A + A_s;
T_a_rep = ones(N,1)*V_s' + V_s*ones(1,N) - 2*T_a_att;

T_a = A.*T_a_att./(V_s*V_s');
T_d = A.*T_a_rep./(V_s*V_s');

W = alpha*T_a - delta*T_d;
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
m = compute_modularity(A,C_out);
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
if(ind_max > 1 && ind_max < N && br_it < N_it)% Is there a new cluster border?
    A_1 = A(1:ind_max,1:ind_max); Deg_1 = Deg(1:ind_max); N_1 = ind_max;
    A_2 = A(ind_max+1:N,ind_max+1:N); Deg_2 = Deg(ind_max+1:N); N_2 = N - ind_max;
    C = [ident_clusters_v2_N(A_1,N_1,Deg_1,L_2,br_it+1,N_it),ind_max,...
         ind_max + ident_clusters_v2_N(A_2,N_2,Deg_2,L_2,br_it+1,N_it)];
else
    C = [];
end;end

function [M,C_vec,N_cl] = optimise_clusters_v1_N(M,C_vec,N_cl,c)
Test = diag(M,1);
while(length(Test) > c-1)
%     disp('Merger')
    % Find position of the largest positive off diagonal element of M
    ind = find(Test == max(Test));
    for br = ind + 1 : N_cl % Update the cluster membership vector
        C_vec(C_vec == br) = br - 1;
    end
    N_cl = N_cl - 1;        % Update number of clusters
    % Update the cluster membership matrix M
    M(:,ind) = M(:,ind) + M(:,ind+1); M(:,ind+1) = [];
    M(ind,:) = M(ind,:) + M(ind+1,:); M(ind+1,:) = [];
    Test = diag(M,1);       % Update the second diagonal vector
end;end

function [M,m,C_out,N_clusters] = ident_communities_v3(A)
% Estimate communities in the network
N = size(A,1);
Deg = A*ones(N,1);
% Step 1: Find cluster borders
L_2 = nnz(A);
C = ident_clusters_v2(A,N,Deg,L_2,0);
% C = ident_clusters_v2(A,N,[],0.18);
if(isempty(C))
    M = 0;m = 0;C_out=ones(N,1);N_clusters = 1;return
end
A_mat = A - (1/nnz(A)).*(Deg*Deg');
%% Step 2: Compute the modularity matrix M
[M,C_out,N_clusters] = compute_modularity_matrix_v1(C,A_mat,A,N);
%% Step 3: Optimise partition by maximising modularity
[M,C_out,N_clusters] = optimise_clusters_v1(M,C_out,N_clusters);
% m = trace(M)
m = compute_modularity(A,C_out);
end

function C = ident_clusters_v2(A,N,Deg,L_2,trsh)
Mod_1 = zeros(N,1);         Mod_2 = zeros(N,1);
Mod_1(1) = -Deg(1)^2/(L_2);   Mod_2(N) = -Deg(N)^2/(L_2);
A_mat = A - (Deg*Deg')./L_2;
for br = 2 : N
    Mod_1(br) = Mod_1(br-1) + 2*sum(A_mat(1:br-1,br)) + A_mat(br,br);
    Mod_2(N-br+1) = Mod_2(N-br+2) + 2*sum(A_mat(N-br+2:N,N-br+1)) + A_mat(N-br+1,N-br+1);
end
[m_max,ind_max] = max(Mod_1+Mod_2);       % Determine cluster border
m_max = m_max(1); ind_max = ind_max(1);
if(m_max >= trsh && ind_max > 1 && ind_max < N)% Is there a new cluster border?
    A_1 = A(1:ind_max,1:ind_max); Deg_1 = Deg(1:ind_max); N_1 = ind_max;
    A_2 = A(ind_max+1:N,ind_max+1:N); Deg_2 = Deg(ind_max+1:N); N_2 = N - ind_max;
    C = [ident_clusters_v2(A_1,N_1,Deg_1,L_2,Mod_1(ind_max)),ind_max,...
         ind_max + ident_clusters_v2(A_2,N_2,Deg_2,L_2,Mod_2(ind_max))];
else
    C = [];
end;end

function [M,C_vec,N_cl] = optimise_clusters_v1(M,C_vec,N_cl)
Test = diag(M,1);
while((max(Test) > 0))
%     disp('Merger')
    % Find position of the largest positive off diagonal element of M
    ind = find(Test == max(Test));
    for br = ind + 1 : N_cl % Update the cluster membership vector
        C_vec(C_vec == br) = br - 1;
    end
    N_cl = N_cl - 1;        % Update number of clusters
    % Update the cluster membership matrix M
    M(:,ind) = M(:,ind) + M(:,ind+1); M(:,ind+1) = [];
    M(ind,:) = M(ind,:) + M(ind+1,:); M(ind+1,:) = [];
    Test = diag(M,1);       % Update the second diagonal vector
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

function N_cl = Non_back_tracking_LCP(A,N,alpha)
Deg = A*ones(N,1);
H = [eye(N) + alpha*(A.*A^2+A) - diag(alpha*(A.*A^2+A)*ones(N,1)) - (eye(N) - diag(Deg)), (eye(N) - diag(Deg)); eye(N), zeros(N)];
Eigs = eigs(H,2*N);
eig_max = maxk(real(Eigs),2);
N_cl = nnz(intersect(find(real(Eigs) > sqrt(eig_max(1))),find(imag(Eigs) == 0)));
end

% Iplementation : Antoine Scherrer
% antoine.scherrer@ens-lyon.fr
% Apply clustering after :
% "Fast unfolding of community hierarchies in large networks"
% Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
% Etienne Lefebvre
% http://arxiv.org/abs/0803.0476
%
% NON ORIENTED VERSION USING SYMETRIC MATRIX A = M + M^t INSTEAD OF 
% (POSSIBLY NON SYMETRIC) INPUT MATRIX M
%
% FULL MATLAB VERSION (SLOWER)
%
% Inputs : 
% M : weight matrix (the matrix is symetrized with
% the sum of weights in both directions)
% s : 1 = Recursive computation
%   : 0 = Just one level computation
% self : 1 = Use self weights
%        0 = Do not use self weights
% debug   : 1 = outputs some debug messages
% verbose : 1 = outputs some messages
%
% Output :
% COMTY, structure with the following information
% for each level i :
%   COMTY.COM{i} : vector of community IDs (sorted by community sizes)
%   COMTY.SIZE{i} : vector of community sizes
%   COMTY.MOD(i) : modularity of clustering
%   COMTY.Niter(i) : Number of iteration before convergence
%
function [COMTY ending] = cluster_jl(M,s,self,debug,verbose)

if nargin < 1
  error('not enough argument');
end

if nargin < 2
  s = 1;
end

if nargin < 3
  self = 1;
end

if nargin < 4
  debug = 0;
end

if nargin < 5
  verbose = 0;
end

S = size(M);
N = S(1);

ddebug = 0;
ending = 0;

% Symetrize matrix taking the sum of weights
M = M + M';
if (self == 0)
  M((N+1).*[0:N-1]+1) = 0;
end
M2 = M;
M2((N+1).*[0:N-1]+1) = 0;

m = sum(sum(M));
Niter = 1;

if m==0 | N == 1
  fprintf('No more possible decomposition\n');
  ending = 1;
  COMTY = 0;
  return;
end

% Main loop
K = sum(M); % Sum of wieght incident to node i
SumTot = sum(M);
SumIn = diag(M); % Sum of weight inside community i
COM = 1:S(1); % Community of node i
for k=1:N
  Neighbor{k} = find(M2(k,:));
end

sCost = 10;
gain = 1;
while (gain == 1)
  Cost = zeros(1,N);
  gain = 0;
  for i=1:N
    Ci = COM(i);
    NB = Neighbor{i};
    G = zeros(1,N); % Gain vector
    best_increase = -1;
    Cnew = Ci;
    COM(i) = -1;
    SumTot(Ci) = SumTot(Ci) - K(i);
    CNj1 = find(COM==Ci);
    SumIn(Ci) = SumIn(Ci) - 2*sum(M(i,CNj1)) - M(i,i);
    for j=1:length(NB)
      Cj = COM(NB(j));
      if (G(Cj) == 0)
        CNj = find(COM==Cj);
        Ki_in = 2*sum(M(i,CNj));
        G(Cj) = Ki_in/m - 2*K(i)*SumTot(Cj)/(m*m);
        if (ddebug)
          fprintf('Gaim for comm %d => %g\n',Cj-1,G(Cj));
        end
        if G(Cj) > best_increase;
          best_increase = G(Cj);
          Cnew_t = Cj;
        end
      end
    end
    if best_increase > 0
      Cnew = Cnew_t;
      if (debug)
        fprintf('Move %d => %d\n',i-1,Cnew-1);
      end
      Cost(i) = best_increase;
    end
    Ck = find(COM==Cnew);
    SumIn(Cnew) = SumIn(Cnew) + 2*sum(M(i,Ck));
    SumTot(Cnew) = SumTot(Cnew) + K(i);
    COM(i) = Cnew;
    if (Cnew ~= Ci)
      gain = 1;
    end
    
  end
  sCost = sum(Cost);
  [C2 S2] = reindex_com(COM);
  Nco = length(unique(COM));
  Nco2 = length(S2(S2>1));
  mod = compute_modularity_jl(COM,M);
  if (debug)
    fprintf('It %d - Mod=%f %d com (%d non isolated)\n',Niter,mod,Nco,Nco2);
  end
  Niter = Niter + 1;
end

Niter = Niter - 1;
[COM COMSIZE] = reindex_com(COM);
COMTY.COM{1} = COM;
COMTY.SIZE{1} = COMSIZE;
COMTY.MOD(1) = compute_modularity_jl(COM,M);
COMTY.Niter(1) = Niter;

% Perform part 2
if (s == 1)
  
  Mnew = M;
  Mold = Mnew;
  COMcur = COM;
  COMfull = COM;
  k = 2;

  if (debug)
    Nco2 = length(COMSIZE(COMSIZE>1));
    fprintf('Pass number 1 - %d com (%d iterations)\n',Nco2,Niter);
  end
  while 1
    Mold = Mnew;
    S2 = size(Mold);
    Nnode = S2(1);
    
    COMu = unique(COMcur);
    Ncom = length(COMu);
    ind_com = zeros(Ncom,Nnode);
    ind_com_full = zeros(Ncom,N);

    for p=1:Ncom
      ind = find(COMcur==p);
      ind_com(p,1:length(ind)) = ind;
    end
    for p=1:Ncom
      ind = find(COMfull==p);
      ind_com_full(p,1:length(ind)) = ind;
    end
    
    Mnew = zeros(Ncom,Ncom);
    for m=1:Ncom
      for n=m:Ncom
        ind1 = ind_com(m,:);
        ind2 = ind_com(n,:);
        Mnew(m,n) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
        Mnew(n,m) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
      end
    end
    
    [COMt e] = cluster_jl(Mnew,0,self,debug,verbose);
    if (e ~= 1)
      COMfull = zeros(1,N);
      COMcur = COMt.COM{1};
      for p=1:Ncom
        ind1 = ind_com_full(p,:);
        COMfull(ind1(ind1>0)) = COMcur(p);
      end
      [COMfull COMSIZE] = reindex_com(COMfull);
      COMTY.COM{k} = COMfull;
      COMTY.SIZE{k} = COMSIZE;
      COMTY.MOD(k) = compute_modularity_jl(COMfull,M);
      COMTY.Niter(k) = COMt.Niter;
      Nco2 = length(COMSIZE(COMSIZE>1));
      if (debug)
        fprintf('Pass number %d - %d com\n',k,Nco2);
      end
      Ind = (COMfull == COMTY.COM{k-1});
      if (sum(Ind) == length(Ind))
        if (debug)
          fprintf('Identical segmentation => End\n');
        end
        return;
      end
    else
      if (debug)
        fprintf('Empty matrix => End\n');
      end
      return;
    end
    k = k + 1;
  end
end

end

% Re-index community IDs
function [C Ss] = reindex_com(COMold)

C = zeros(1,length(COMold));
COMu = unique(COMold);
S = zeros(1,length(COMu));
for l=1:length(COMu)
    S(l) = length(COMold(COMold==COMu(l)));
end
[Ss INDs] = sort(S,'descend');

for l=1:length(COMu)
    C(COMold==COMu(INDs(l))) = l;
end

end

%Compute modulartiy
function MOD = compute_modularity_jl(C,Mat)

m = sum(sum(Mat));
MOD = 0;
COMu = unique(C);
for j=1:length(COMu)
    Cj = find(C==COMu(j));
    Ec = sum(sum(Mat(Cj,Cj)));
    Et = sum(sum(Mat(Cj,:)));
    if Et>0
        MOD = MOD + Ec/m-(Et/m)^2;
    end
end

end

% This function implements the algorithm of Newton
function [N_cl,m] = Newman_clustering(A)
% clc
N = size(A,1);
Deg = A*ones(N,1);
M = A - (Deg*Deg')./nnz(A);
C = ones(N,1);
C_out = Newman_ident_cl(C,M,0);
% [val,ind] = sort(C_out);
% plot(sort(C_out))
% figure
% spy(A(ind,ind))
% pause()
N_cl = length(unique(C_out));
m = compute_modularity(A,C_out);
end

function C_out = Newman_ident_cl(C,M,tr)
C_out = C;
[X,L] = eig(M);
[val,ind] = max(diag(L));
ind_pos = find(X(:,ind) > 0);
ind_neg = find(X(:,ind) < 0);

if(val <= 0 || isempty(ind_pos) || isempty(ind_neg))
    return
else   
    C_out(ind_pos) = C(ind_pos).*rand();
    C_out(ind_neg) = C(ind_neg).*rand();

    M_pos = M(ind_pos,ind_pos);
    M_neg = M(ind_neg,ind_neg);
    
    m_pos = sum(sum(M_pos));
    m_neg = sum(sum(M_neg));

    if(m_pos + m_neg <= tr)
        return
    else
        M_pos = M_pos - diag(M_pos*ones(length(ind_pos),1));
        M_neg = M_neg - diag(M_neg*ones(length(ind_neg),1));
   
%         C_out(ind_pos) = Newman_ident_cl(C_out(ind_pos),M_pos,m_pos);
%         C_out(ind_neg) = Newman_ident_cl(C_out(ind_neg),M_neg,m_neg);
        C_out(ind_pos) = Newman_ident_cl(C_out(ind_pos),M_pos,0);
        C_out(ind_neg) = Newman_ident_cl(C_out(ind_neg),M_neg,0);
    end  
end
end

function N_cl = Non_back_tracking(links_dir,A,N)
% L_2 = size(links_dir,1);
% B = zeros(L_2);
% for i = 1 : L_2
%     orig = links_dir(i,1);
%     dest = links_dir(i,2);
%     ind_1 = find(links_dir(:,1) == dest);
%     ind_2 = find(links_dir(ind_1,2)~=orig);
%     B(i,ind_1(ind_2)) = 1;
% end
% Eigs = eig(B);
% % plot(real(Eigs),imag(Eigs),'.')
% ind = find(imag(Eigs) == 0);
% N_cl = sum(Eigs(ind) > sqrt(max(Eigs)));
B_alt = [A, eye(N) - diag(A*ones(N,1));eye(N) zeros(N)];
Eigs = eigs(B_alt,2*N);
% plot(real(Eigs),imag(Eigs),'.')
ind = find(imag(Eigs) == 0);
N_cl = sum(Eigs(ind) > sqrt(max(Eigs)));
end
