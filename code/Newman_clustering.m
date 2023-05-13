function [C_out,N_cl,m] = Newman_clustering(A)
N = size(A,1);
Deg = A*ones(N,1);
M = A - (Deg*Deg')./nnz(A);
C = ones(N,1);
C_out = Newman_ident_cl(C,M,0);
m = compute_modularity(A,C_out);
uniComm = unique(C_out);
for i=1:length(uniComm)
    C_out(C_out==uniComm(i)) = max(uniComm) +i;
end
uniComm = unique(C_out);
for i=1:length(uniComm)
    C_out(C_out==uniComm(i)) = i;
end
N_cl = max(C_out);
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

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end