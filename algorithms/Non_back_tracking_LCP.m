function [Comm,N_cl,Q] = Non_back_tracking_LCP(A,N)                                   % LCP-based non back tracking approach for estimating number of clusters (Section 4.4)
Deg = A*ones(N,1);                                                                          % Degree vector
gamma = fit_powerlaw(Deg); %fit power-law exponent gamma
if gamma <= 1 %determine if the network conform power-law
    alpha = 0.95;
elseif gamma <= 3.5
    alpha = 0.3*gamma - 0.1; 
else
    alpha = 0.95;
end

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

function gamma = fit_powerlaw(degree)
prob = tabulate(degree);
prob(:,2) = [];
prob_0 = prob(:,2) == 0;
prob(prob_0,:) = [];
prob = log(prob);
gamma = -polyfit(prob(:,1),prob(:,2),1)*[1 0]';
end
