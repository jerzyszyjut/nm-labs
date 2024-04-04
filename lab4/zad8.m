a = 1;
b = 60000;
ytolerance = 1e-12;
max_iterations = 100;

[n_bisection, ~, ~, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @(N) estimate_execution_time(N));
[n_secant, ~, ~, xtab_secant, xdif_secant] = secant_method(a, b, max_iterations, ytolerance, @(N) estimate_execution_time(N));

figure;
subplot(2, 1, 1);
plot(xtab_bisection, 'r');
hold on;
plot(xtab_secant, 'b');
title('Przybliżone wartości N');
xlabel('Iteracja');
ylabel('t');
legend('Metoda bisekcji', 'Metoda siecznych');

subplot(2, 1, 2);
semilogy(xdif_bisection, 'r');
hold on;
semilogy(xdif_secant, 'b');
title('Różnica pomiędzy kolejnymi przybliżeniami N');
xlabel('Iteracja');
ylabel('Różnica');
legend('Metoda bisekcji', 'Metoda siecznych');

saveas(gcf, 'zadanie_8.png');

function time_delta = estimate_execution_time(N)
    % time_delta - różnica pomiędzy estymowanym czasem wykonania algorytmu dla zadanej wartości N a zadanym czasem M
    % N - liczba parametrów wejściowych
    if N <= 0
        error('N must be greater than 0');
    end
    M = 5000; % [s]
    
    time_delta = abs((N^(16/11) + N^((pi^2)/8))/1000) - M;
end

function [xsolution,ysolution,iterations,xtab,xdif] = bisection_method(a,b,max_iterations,ytolerance,fun)
    % a - lewa granica przedziału poszukiwań miejsca zerowego
    % b - prawa granica przedziału poszukiwań miejsca zerowego
    % max_iterations - maksymalna liczba iteracji działania metody bisekcji
    % ytolerance - wartość abs(fun(xsolution)) powinna być mniejsza niż ytolerance
    % fun - nazwa funkcji, której miejsce zerowe będzie wyznaczane
    %
    % xsolution - obliczone miejsce zerowe
    % ysolution - wartość fun(xsolution)
    % iterations - liczba iteracji wykonana w celu wyznaczenia xsolution
    % xtab - wektor z kolejnymi kandydatami na miejsce zerowe, począwszy od xtab(1)= (a+b)/2
    % xdiff - wektor wartości bezwzględnych z różnic pomiędzy i-tym oraz (i+1)-ym elementem wektora xtab; xdiff(1) = abs(xtab(2)-xtab(1));
    iterations = 0;
    xtab = [];
    xdif = [];
    error = 1;

    while iterations < max_iterations && error > ytolerance
        c = (a+b)/2;
        xtab = [xtab; c];
        if fun(a)*fun(c) < 0
            b = c;
        else
            a = c;
        end
        iterations = iterations + 1;
        if iterations > 1
            xdif = [xdif; abs(xtab(end)-xtab(end-1))];
            error = abs(fun(c));
        end
    end
    
    xsolution = c;
    ysolution = fun(c);
end

function [xsolution,ysolution,iterations,xtab,xdif] = secant_method(a,b,max_iterations,ytolerance,fun)
    % a - lewa granica przedziału poszukiwań miejsca zerowego (x0=a)
    % b - prawa granica przedziału poszukiwań miejsca zerowego (x1=b)
    % max_iterations - maksymalna liczba iteracji działania metody siecznych
    % ytolerance - wartość abs(fun(xsolution)) powinna być mniejsza niż ytolerance
    % fun - nazwa funkcji, której miejsce zerowe będzie wyznaczane
    %
    % xsolution - obliczone miejsce zerowe
    % ysolution - wartość fun(xsolution)
    % iterations - liczba iteracji wykonana w celu wyznaczenia xsolution
    % xtab - wektor z kolejnymi kandydatami na miejsce zerowe, począwszy od x2
    % xdiff - wektor wartości bezwzględnych z różnic pomiędzy i-tym oraz (i+1)-ym elementem wektora xtab; xdiff(1) = abs(xtab(2)-xtab(1));
    iterations = 0;
    xtab = [a;b];
    xdif = [];
    error = 1;

    while iterations < max_iterations && error > ytolerance
        x_k = xtab(end);
        x_k_1 = xtab(end-1);
        c = x_k - (fun(x_k)*(x_k - x_k_1)/(fun(x_k) - fun(x_k_1)));
        xtab = [xtab; c];

        iterations = iterations + 1;
        if iterations > 1
            xdif = [xdif; abs(xtab(end)-xtab(end-1))];
            error = abs(fun(c));
        end
    end
    xtab = xtab(3:end);

    xsolution = c;
    ysolution = fun(c);
end
