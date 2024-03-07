function main
    [numer_indeksu, Edges, I, B, A, b, M, r] = page_rank();
    bar(r);
    title('PageRank');
    xlabel('Numer strony');
    ylabel('Wartość PageRank');
    print -dpng 'zadanie7.png';
end

function [numer_indeksu, Edges, I, B, A, b, M, r] = page_rank()
    N = 8;
    d = 0.85;
    numer_indeksu = 193064;
    L1 = mod(floor(numer_indeksu/10), 10);
    L2 = mod(floor(numer_indeksu/100), 10);
    do = mod(L1, 7) + 1;
    z = mod(L2, 7) + 1;
    Edges = [   1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, z;
                4, 6, 3, 4, 5, 5, 6, 7, 5, 6, 4, 6, 4, 7, 6, do, 8];
    
    I = speye(N);
    
    B = sparse(Edges(2, :), Edges(1, :), 1, N, N);

    A = spdiags(1./sum(B)', 0, N, N);
    proportion = ((1-d)/N);
    b = (proportion * ones(N, 1));
    M = I-(d*B*A);
    r = M \ b;
end