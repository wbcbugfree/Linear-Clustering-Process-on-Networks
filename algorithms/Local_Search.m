function [Comm,c,Q]=Local_Search(A)
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
LL_D = D(C);
max_LL_D = max(LL_D);
index_LL_D = [];
for i = 1:length(LL_D)
    if LL_D(i)==max_LL_D
        index_LL_D = [index_LL_D i];
    end
end
M = C(index_LL_D); %local leaders with max degree

%links towards local leader with larger degree and shortest path
[~,index_sort] = sort(LL_D);
C_sort = C(index_sort);
for i = 1:length(C)-length(M)
    Distance = [];
    for j = i+1:length(C)
        [~,dis] = shortestpath(G,C_sort(i),C_sort(j));
        Distance = [Distance dis];
    end
    Dist_flip = fliplr(Distance);
    [min_Dis,ind_dis_flip] = min(Dist_flip);
    l(C_sort(i)) = min_Dis;
    ind_dis = length(Distance)-ind_dis_flip+1;
    diG_u = [diG_u; [C_sort(i), C_sort(ind_dis+i)]];
end

for i = 1:length(M)
    l(M(i)) = max(l);
end

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

for i = 1:length(C)
    delta(i) = ((l_star(C(i))-min(l_star))/(max(l_star)-min(l_star)))*((D_star(C(i))-min(D_star))/(max(D_star)-min(D_star)));
end

%assign community labels
Comm = [];
for i = 1:N
    if ismember(i,C)
        Comm(i) = i;
    else
        Comm(i) = diG_u(diG_u(:,1)==i,2);
    end
end

while(~isequal(sort(unique(Comm)),sort(C)))
    for i = 1:N
        if ~ismember(Comm(i),C)
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
