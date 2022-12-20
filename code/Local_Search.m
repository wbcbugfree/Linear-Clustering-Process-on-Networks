function [delta]=Local_Search(G)

D = degree(G);
A = adjacency(G);
N = size(A,1);
diG=[0,0];
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
    for j = 1:length(max_Nb) %build digraph
        if D(max_Nb(j))>=D(i)&&max(ismember(diG,[max_Nb(j),i],'rows'))==0
            diG = [diG; [i, max_Nb(j)]];
            l(i)=1;
        end
    end
end
diG(1,:) = [];

%remove multiple out-going links
uni = unique(diG(:,1));
hist = histc(diG(:,1),uni);
repeat = find(hist>=2);
while(~isempty(repeat))
    diG(repeat(1),:) = [];
    repeat(1) = [];    
    uni = unique(diG(:,1));
    hist = histc(diG(:,1),uni);
    repeat = find(hist>=2);
end

edgetable = table([diG(:,1) diG(:,2)]);
edgetable.Properties.VariableNames(1) = "EndNodes";
DG = digraph(edgetable);
OD = outdegree(DG);
C = find(OD==0); %local leaders
LL_D = D(C);
max_LL_D = max(LL_D);
index_LL_D = [];
for i = 1:length(LL_D)
    if LL_D(i)==max_LL_D
    index_LL_D = [index_LL_D i];
    end
end
M = C(index_LL_D); %local leaders with max degree

%links between local leaders
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
    diG = [diG; [C_sort(i), C_sort(ind_dis+i)]];
end
edgetable = table([diG(:,1) diG(:,2)]);
edgetable.Properties.VariableNames(1) = "EndNodes";
DG = digraph(edgetable);
plot(DG)

for i = 1:length(M)
    l(M(i)) = max(l);
end

D_star = D;
D_sorted = sort(D);
D_uni = unique(D_sorted);
for i = 1:length(D_star)
    for j = 1:length(D_uni)
        if D_star(i)==D_uni(j)
            D_star(i) = j;
        end
    end
end

l_star = l.^2;

for i = 1:length(C)
    delta(i) = ((l_star(C(i))-min(l_star))/(max(l_star)-min(l_star)))*((D_star(C(i))-min(D_star))/(max(D_star)-min(D_star)));
end
delta = sort(delta, 'descend');

end
