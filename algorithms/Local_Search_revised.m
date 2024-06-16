function [Comm,c,Q]=Local_Search_revised(A)
G = graph(A);
D = degree(G);
N = size(A,1);
diG = [0,0];
l = zeros(1,N);

for i = 1:N
    Nb = neighbors(G,i); %find neighbors
    D_Nb = D(Nb); %degree of neighbors
    max_D_Nb = max(D_Nb); %max degree of neighbors
    index_D_Nb = []; %index of max
    for j = 1:length(D_Nb) %return the index, maybe more than 1
        if D_Nb(j)==max_D_Nb
            index_D_Nb = [index_D_Nb j];
        end
    end
    max_Nb = Nb(index_D_Nb); %neighbors with max degree
    for k = 1:length(max_Nb) %build digraph
        if D(max_Nb(k))>=D(i)&&max(ismember(diG,[max_Nb(k),i],'rows'))==0
            diG = [diG; [i, max_Nb(k)]];
            l(i)=1;
        end
    end
end
diG(1,:) = [];

%remove multiple out-going links
uni = unique(diG(:,1));
diG_u = [0,0];
for i = 1:length(uni)
    for j = 1:length(diG(:,1))
        if uni(i) == diG(j,1)
            diG_u(i,1) = uni(i);
            diG_u(i,2) = diG(j,2);
            break
        end
    end
end

C = setdiff(linspace(1,N,N),diG_u(:,1)); %local leaders
D_C = D(C);
[~,index_C] = sort(D_C, 'descend');
C_des = C(index_C);
C_new = [];
if length(C) == 1
    C_new = C;
    l(C) = max(l);
else
    for i = 1:length(C)%links between local leaders
        dest_set = setdiff(C_des, C(i), 'stable');
        Distance = [];
        for j = 1:length(dest_set)
            [~, dis] = shortestpath(G,C(i),dest_set(j));
            Distance = [Distance dis];
        end
        [min_Dis,ind_dis] = min(Distance);
        l(C(i)) = min_Dis;
        diG_u = [diG_u; [C(i), dest_set(ind_dis)]];
        C_new = [C_new dest_set(ind_dis)];
    end
end

C_new = unique(C_new);

D_star = D;
D_sorted = sort(D);
D_uni = unique(D_sorted);
for i = 1:length(D_star)
    for j = 1:length(D_uni)
        if D_star(i) == D_uni(j)
            D_star(i) = j;
        end
    end
end

l_star = l.^2;

for i = 1:length(C_new)
    delta(i) = ((l_star(C_new(i))-min(l_star))/(max(l_star)-min(l_star)))*((D_star(C_new(i))-min(D_star))/(max(D_star)-min(D_star)));
end

%assign community labels
Comm = [];
for i = 1:N
    if ismember(i,C_new)
        Comm(i) = i;
    else
        Comm(i) = diG_u(diG_u(:,1)==i,2);
    end
end

while(~isequal(sort(unique(Comm)),sort(C_new)))
    for i = 1:N
        if ~ismember(Comm(i),C_new)
            Comm(i) = diG_u(diG_u(:,1)==Comm(i),2);
        end
    end
end

uniComm = unique(Comm);
for i=1:length(uniComm)
    Comm(Comm==uniComm(i)) = max(uniComm) + i;
end
uniComm = unique(Comm);
for i=1:length(uniComm)
    Comm(Comm==uniComm(i)) = i;
end
c = max(Comm);
Q = compute_modularity(A,Comm');
end

function Q = compute_modularity(A,C)
N = size(A,1);
Deg = A*ones(N,1);
Q = 1/nnz(A).*sum(sum((A - (1./nnz(A).*Deg*Deg')).*((ones(N,1)*abs(C)' - C*ones(1,N)) == 0)));
end
