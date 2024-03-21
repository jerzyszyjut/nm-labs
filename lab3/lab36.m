function main()
    load('filtr_dielektryczny.mat');
    iterations = 500;
    [err_norm_J, err_norms_J, time_J] = solve_Jacobi(A, b, iterations);
    time_J;
    err_norm_J;
    [err_norm_GS, err_norms_GS, time_GS] = solve_Gauss_Seidel(A, b, iterations);
    time_GS;
    err_norm_GS;
    [~, time] = solve_direct(A, b);
    time;

    figure;
    plot(1:length(err_norms_GS), err_norms_GS, 'b');
    title('Norma błędu rezydualnego w zależności od liczby iteracji dla Gauss-Seidela');
    xlabel('Liczba iteracji');
    ylabel('Norma błędu rezydualnego');
    legend('Gauss-Seidel');
    hold off;
    saveas(gcf, 'zasanie6GSplot.png');

    figure;
    semilogy(1:length(err_norms_J), err_norms_J, 'r');
    title('Norma błędu rezydualnego w zależności od liczby iteracji dla Jacobiego');
    xlabel('Liczba iteracji');
    ylabel('Norma błędu rezydualnego');
    legend('Jacobi');
    hold off;
    saveas(gcf, 'zasanie6Jplot.png');

end

function [err_norm, err_norms, time] = solve_Gauss_Seidel(A, b, iterations_limit)
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b)
    % err_norms - wektor norm błędu rezydualnego rozwiązania x w kolejnych iteracjach
    % time - czas wyznaczenia rozwiązania x
        
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    N = length(A);
    x = ones(N,1);
    err_norms = [];
    iterations = 0;
    
    M = -(D+L)\U;
    bm = (D+L)\b;
    
    tic;
    err_norm = norm(A*x-b);
    while err_norm > 1e-5 && iterations < iterations_limit
        x = M*x + bm;
        err_norm = norm(A*x-b);
        err_norms = [err_norms, err_norm];
        iterations = iterations + 1;
    end
    time = toc;
end

function [err_norm, err_norms, time] = solve_Jacobi(A, b, iterations_limit)
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b);
    % err_norms - wektor norm błędu rezydualnego rozwiązania x w kolejnych iteracjach
    % time - czas wyznaczenia rozwiązania x
    % index_number - Twój numer indeksu
    
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    N = length(A);
    x = ones(N,1);
    M = -D\(L+U);
    bm = D\b;
    iterations = 0;
    err_norms = [];

    tic;
    err_norm = norm(A*x-b);
    while norm(A*x-b) > 1e-5 && iterations < iterations_limit
        x = M*x + bm;
        err_norm = norm(A*x-b);
        err_norms = [err_norms, err_norm];
        iterations = iterations + 1;
    end
    time = toc;
end

function [time,err_norm] = solve_direct(A,b)
    % time - czas wyznaczenia rozwiązania x
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b);

    tic;
    x = A\b;
    time = toc;
    err_norm = norm(A*x-b);
end
