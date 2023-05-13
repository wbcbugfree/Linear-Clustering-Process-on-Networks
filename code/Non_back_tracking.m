function [Comm,N_cl,Q] = Non_back_tracking(A,N)
B_star = [A, eye(N) - diag(A*ones(N,1));eye(N) zeros(N)];       % Compute the 2N x 2N matrix B_star in eq. (27), page 14
Eigs = eigs(B_star,2*N);                                        % Compute the eigenvalues of B_star
ind = imag(Eigs) == 0;                                          % Identify real eigenvalues
N_cl = sum(Eigs(ind) > sqrt(max(Eigs)));                        % Determine the number of real eigenvalues larger that sqrt(lambda_1)    
[U,~] = eigs(B_star,N_cl,'la');
U = real(U);
Comm = kmeans(U(1:N,:),N_cl);
Q = compute_modularity(A,Comm);
end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end
