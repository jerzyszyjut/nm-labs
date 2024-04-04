a = 1;
b = 50;
ytolerance = 1e-12;
max_iterations = 1000;
omega = 75;

[omega_bisection, ~, ~, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @(omega) impedance_magnitude(omega));
[omega_secant, ~, ~, xtab_secant, xdif_secant] = secant_method(a, b, max_iterations, ytolerance, @(omega) impedance_magnitude(omega));

figure;
subplot(2, 1, 1);
plot(xtab_bisection, 'r');
hold on;
plot(xtab_secant, 'b');
title('Przybliżone wartości omegi');
xlabel('Iteracja');
ylabel('Omega');
legend('Metoda bisekcji', 'Metoda siecznych');

subplot(2, 1, 2);
semilogy(xdif_bisection, 'r');
hold on;
semilogy(xdif_secant, 'b');
title('Różnica pomiędzy kolejnymi przybliżeniami omegi');
xlabel('Iteracja');
ylabel('Różnica');
legend('Metoda bisekcji', 'Metoda siecznych');

saveas(gcf, 'zadanie_4.png');

function impedance_delta = impedance_magnitude(omega)
    if omega <= 0
        error('Omega must be greater than 0');
    end
    
    R = 525;
    C = 7e-5;
    L = 3;
    M = 75; % docelowa wartość modułu impedancji
    
    Z = abs(1/sqrt(1/(R^2) + (omega * C - 1/(omega * L)) ^ 2));
    impedance_delta = Z - M;
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