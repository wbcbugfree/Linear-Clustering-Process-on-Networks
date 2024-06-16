function [Comm,c,Q] = Eigengap_B(A)
N = size(A,1);
Deg = A*ones(N,1);
L_2 = sum(sum(A));
M = A - Deg*Deg'./L_2;
L = eigs(M,N);                              % Compute eigenvalues of Q
Eigs = sort(L,'ascend');                    % Sort the eigenvalues of Q
Eig_gap = flipud(diff((Eigs)));
[~,c] = max(Eig_gap(1:ceil(0.5*N)));        % Determine the largest gap
c = c + 1;
[U,~] = eigs(M,c,'la');
Comm = kmeans(U,c);
Q = compute_modularity(A,Comm);
end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end
