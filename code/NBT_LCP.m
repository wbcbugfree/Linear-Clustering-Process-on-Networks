function [Comm,N_cl,Q] = NBT_LCP(A,N,alpha)                                   % LCP-based non back tracking approach for estimating number of clusters (Section 4.4)
Deg = A*ones(N,1);                                                                          % Degree vector
H = [eye(N) + alpha*(A.*A^2+A) - diag(alpha*(A.*A^2+A)*ones(N,1)) - (eye(N) - diag(Deg)), (eye(N) - diag(Deg)); eye(N), zeros(N)];  % The 2Nx2N matrix in eq. (28), page 15
Lambda_H = eigs(H,2*N);                                                                     % Compute the eigenvalues of H
eig_max = maxk(real(Lambda_H),2);                                                           % Determine the largest eigenvalue
N_cl = nnz(intersect(find(real(Lambda_H) > sqrt(eig_max(1))),find(imag(Lambda_H) == 0)));   % Determine the number of real eigenvalues larger that sqrt(lambda_1)
[U,~] = eigs(H,N_cl,'la');
U = real(U);
Comm = kmeans(U(1:N,:),N_cl);
Q = compute_modularity(A,Comm);
end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end