function main()
    N = 1000:1000:8000;
    plot_problem_5(N);
end

function plot_problem_5(N)
    [time_Jacobi, iterations_Jacobi, time_Gauss_Seidel, iterations_Gauss_Seidel] = benchmark(N);
    figure;
    subplot(2,1,1);
    plot(N,time_Jacobi);
    hold on;
    plot(N,time_Gauss_Seidel);
    xlabel('Rozmiar macierzy');
    ylabel('Czas wyznaczenia rozwiązania');
    title('Czas wyznaczenia rozwiązania w zależności od rozmiaru macierzy');
    legend('Jacobi','Gauss-Seidel', 'Location', 'eastoutside');
    hold off;
    subplot(2,1,2);
    bar(N,iterations_Jacobi);
    hold on;
    bar(N,iterations_Gauss_Seidel);
    xlabel('Rozmiar macierzy');
    ylabel('Liczba iteracji');
    title('Liczba iteracji w zależności od rozmiaru macierzy');
    legend('Jacobi','Gauss-Seidel', 'Location', 'eastoutside');
    hold off;
    saveas(gcf, 'zadanie5.png');
end

function [time_Jacobi, iterations_Jacobi, time_Gauss_Seidel, iterations_Gauss_Seidel] = benchmark(N)
    % time_Jacobi - czas wyznaczenia rozwiązania równania macierzowego metodą Jacobiego
    % iterations_Jacobi - liczba iteracji wykonana w procesie iteracyjnym metody Jacobiego
    % time_Gauss_Seidel - czas wyznaczenia rozwiązania równania macierzowego metodą Gaussa-Seidla
    % iterations_Gauss_Seidel - liczba iteracji wykonana w procesie iteracyjnym metody Gaussa-Seidla
    time_Jacobi = zeros(1,length(N));
    iterations_Jacobi = zeros(1,length(N));
    time_Gauss_Seidel = zeros(1,length(N));
    iterations_Gauss_Seidel = zeros(1,length(N));

    for i=1:length(N)
        [~,~,~,~,~,~,time_Jacobi(i),iterations_Jacobi(i),~] = solve_Jacobi(N(i));
        [~,~,~,~,~,~,time_Gauss_Seidel(i),iterations_Gauss_Seidel(i),~] = solve_Gauss_Seidel(N(i));
    end
end

function plot_direct(N,vtime_direct)
    for i=1:length(N)
        [~,~,~,time_direct,~,~] = solve_direct(N(i));
        vtime_direct(i) = time_direct;
    end
    figure;
    plot(N,vtime_direct);
    xlabel('Rozmiar macierzy');
    ylabel('Czas wyznaczenia rozwiązania');
    title('Czas wyznaczenia rozwiązania w zależności od rozmiaru macierzy');
    saveas(gcf, 'zadanie2.png');
end

function [A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Gauss_Seidel(N)
    % A - macierz rzadka z równania macierzowego A * x = b
    % b - wektor prawej strony równania macierzowego A * x = b
    % M - macierz pomocnicza opisana w instrukcji do Laboratorium 3 – sprawdź wzór (7) w instrukcji, który definiuje M jako M_{GS}
    % bm - wektor pomocniczy opisany w instrukcji do Laboratorium 3 – sprawdź wzór (7) w instrukcji, który definiuje bm jako b_{mGS}
    % x - rozwiązanie równania macierzowego
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b)
    % time - czas wyznaczenia rozwiązania x
    % iterations - liczba iteracji wykonana w procesie iteracyjnym metody Gaussa-Seidla
    % index_number - Twój numer indeksu
    index_number = 193064;
    L1 = mod(index_number, 10);
        
    [A,b] = generate_matrix(N, L1);

    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    x = ones(N,1);
    iterations = 0;
    
    M = -(D+L)\U;
    bm = (D+L)\b;
    
    tic;
    while norm(A*x-b) > 1e-5
        x = M*x + bm;
        iterations = iterations + 1;
    end
    time = toc;
    
    err_norm = norm(A*x-b);
end

function [A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Jacobi(N)
    % A - macierz z równania macierzowego A * x = b
    % b - wektor prawej strony równania macierzowego A * x = b
    % x - rozwiązanie równania macierzowego
    % time_direct - czas wyznaczenia rozwiązania x
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b);
    % index_number - Twój numer indeksu
    index_number = 193064;
    L1 = mod(index_number, 10);

    [A,b] = generate_matrix(N, L1);
    

    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);
    x = ones(N,1);
    iterations = 0;
    M = -D\(L+U);
    bm = D\b;

    tic;
    while norm(A*x-b) > 1e-5
        x = M*x + bm;
        iterations = iterations + 1;
    end
    time = toc;

    err_norm = norm(A*x-b);
end

function [A,b,x,time_direct,err_norm,index_number] = solve_direct(N)
    % A - macierz z równania macierzowego A * x = b
    % b - wektor prawej strony równania macierzowego A * x = b
    % x - rozwiązanie równania macierzowego
    % time_direct - czas wyznaczenia rozwiązania x
    % err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b);
    % index_number - Twój numer indeksu
    index_number = 193064;
    L1 = mod(index_number, 10);

    [A,b] = generate_matrix(N, L1);
    tic;
    x = A\b;
    time_direct = toc;
    err_norm = norm(A*x-b);
end

function [A,b] = generate_matrix(N, convergence_factor)
    % A - macierz o rozmiarze NxN
    % b - wektor o rozmiarze Nx1
    % convergense_factor - regulacja elementów diagonalnych macierzy A, które wpływają
    %       na zbieżność algorytmów iteracyjnego rozwiązywania równania macierzowego

    if(convergence_factor<0 || convergence_factor>9)
        error('Wartość convergence_factor powinna być zawarta w przedziale [1,9]');
    end

    seed = 0; % seed - kontrola losowości elementów niezerowych macierzy A
    rng(seed); % ustawienie generatora liczb losowych

    A = rand(N, N);
    A = A - diag(diag(A)); % wyzerowanie głównej diagonalnej

    convergence_factor_2 = 1.2 + convergence_factor/10;
    diag_values = sum(abs(A),2) * convergence_factor_2;
    A = A + diag(diag_values); % nadanie nowych wartości na głównej diagonalnej

    % regulacja normy macierzy
    norm_Frobenius = norm(A,'fro');
    A = A/norm_Frobenius;

    b = rand(N,1);
end

