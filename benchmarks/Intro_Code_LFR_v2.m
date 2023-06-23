clc
clear
close all

[A,C_0,c,m] = LFR_benchmark_v9(500,12,5,2,3,0.25);
N = size(A,1);

[~,pos_0] = sort(C_0);
figure(1)
subplot(2,3,1)
spy(A(pos_0,pos_0))
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
z_LCP = nmi(C_LCP, C_0)
disp('NMI_LOU:')
z_Lou = nmi(C_LOU, C_0)
disp('NMI_Le:')
z_Le = nmi(C_Le, C_0)
disp('NMI_NEW:')
z_New = nmi(C_NEW, C_0)
disp('NMI_nbt:')
z_nbt = nmi(C_nbt, C_0)
% disp('NMI_LCP_n:')
% z_LCP_n = nmi(C_LCP_n, C_0)
% disp('NMI_eig:')
% z_eig = nmi(C_eig, C_0)
% disp('NMI_LS:')
% z_LS = nmi(C_LS, C_0)

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