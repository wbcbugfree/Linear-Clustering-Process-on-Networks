clc
clear
close all

N_clusters = 20;
N = 500;                                                               % Number of nodes
d_av = 7;

% b_in = linspace(8,14,20);   % c = 2 clusters
% b_in = linspace(9,20,20);   % c = 3 clusters
% b_in = linspace(11,25,20);   % c = 4 clusters
% b_in = linspace(21,55,20);   % c = 8 clusters
% b_in = linspace(17,62,20);   % c = 10 clusters
% b_in = linspace(45,125,20);   % c = 20 clusters

b_in = 125;
b_out = (d_av*N_clusters - b_in)./(N_clusters - 1);

[G,G_bidir,links,links_dir,nodes,groups] = generate_sbm(N,N_clusters,b_in,b_out);
A = adjacency(G);
N = size(A,1);
Deg = A*ones(N,1);

ind = find(Deg>0);
A = A(ind,ind);
N = length(ind);
Deg = Deg(ind);

C = groups(ind)';
groups = groups(ind);
[val,pos] = sort(C);
C = C(pos);
A_0 = A(pos,pos);
Q_0 = compute_modularity(A_0,C);

figure(1)
subplot(2,3,1)
spy(A_0)
title('Original LFR graph','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

[C_LCP,c_LCP,m_LCP] = LCP_v2(A);
[~,pos_LCP] = sort(C_LCP);
figure(1)
subplot(2,3,2)
spy(A(pos_LCP,pos_LCP))
title('LCP output','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

[C_LOU,c_LOU,m_LOU] = Louvain_Algorithm_v2(A);
[~,pos_LOU] = sort(C_LOU);
figure(1)
subplot(2,3,3)
spy(A(pos_LOU,pos_LOU))
title('Louvain output','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

[C_Le,c_Le,m_Le] = Leiden_v1(A);
[~,pos_Le] = sort(C_Le);
figure(1)
subplot(2,3,4)
spy(A(pos_Le,pos_Le))
title('Leiden output','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

[C_NEW,c_NEW,m_NEW] = Newman_clustering(A);
[~,pos_NEW] = sort(C_NEW);
figure(1)
subplot(2,3,5)
spy(A(pos_NEW,pos_NEW))
title('Newman output','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

[C_nbt,c_nbt,m_nbt] = Non_back_tracking(A,N);
[~,pos_nbt] = sort(C_nbt);
figure(1)
subplot(2,3,6)
spy(A(pos_nbt,pos_nbt))
title('Non-backtracking output','Interpreter','latex')
set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
axis off

% [C_LCP_n,c_LCP_n,m_LCP_n] = Non_back_tracking_LCP(A,N);
% [~,pos_LCP_n] = sort(C_LCP_n);
% figure(1)
% subplot(3,3,7)
% spy(A(pos_LCP_n,pos_LCP_n))
% title('Non-backtracking LCP output','Interpreter','latex')
% set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
% axis off
% 
% [C_eig,c_eig,m_eig] = Eigengap_B(A);
% [~,pos_eig] = sort(C_eig);
% figure(1)
% subplot(3,3,8)
% spy(A(pos_eig,pos_eig))
% title('Modularity Eigengap output','Interpreter','latex')
% set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
% axis off

% [C_LS,c_LS,m_LS] = Local_Search(A);
% [~,pos_LS] = sort(C_LS);
% figure(1)
% subplot(3,3,9)
% spy(A(pos_LS,pos_LS))
% title('Local Search output','Interpreter','latex')
% set(gca,'Fontsize',12,'TickLabelInterpreter','latex')
% axis off

disp('_____________________NMI_____________________')
disp('NMI_LCP:')
z_LCP = nmi(C_LCP, groups)
disp('NMI_LOU:')
z_Lou = nmi(C_LOU, groups)
disp('NMI_Le:')
z_Le = nmi(C_Le, groups)
disp('NMI_NEW:')
z_New = nmi(C_NEW, groups)
disp('NMI_nbt:')
z_nbt = nmi(C_nbt, groups)
% disp('NMI_LCP_n:')
% z_LCP_n = nmi(C_LCP_n, groups)
% disp('NMI_eig:')
% z_eig = nmi(C_eig, groups)
% disp('NMI_LS:')
% z_LS = nmi(C_LS, groups)

%% Used functions
function z = nmi(x, y)
% Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
% Input:
%   x, y: two integer vector of the same length 
% Ouput:
%   z: normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))
% Written by Mo Chen (sth4nth@gmail.com).
assert(numel(x) == numel(y));
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);
l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));
idx = 1:n;
Mx = sparse(idx,x,1,n,k,n);
My = sparse(idx,y,1,n,k,n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log2(Pxy));
% hacking, to elimative the 0log0 issue
Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));
% entropy of Py and Px
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));
% mutual information
MI = Hx + Hy - Hxy;
% normalized mutual information
z = sqrt((MI/Hx)*(MI/Hy));
z = max(0,z);
end

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
