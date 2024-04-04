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
