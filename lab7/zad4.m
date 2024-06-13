[integration_error, Nt, ft_5, xr, yr, yrmax] = zadanie4();

function [integration_error, Nt, ft_5, xr, yr, yrmax] = zadanie4()
    % Numeryczne całkowanie metodą Monte Carlo.
    %
    %   integration_error - wektor wierszowy. Każdy element integration_error(1,i)
    %       zawiera błąd całkowania obliczony dla liczby losowań równej Nt(1,i).
    %       Zakładając, że obliczona wartość całki dla Nt(1,i) próbek wynosi
    %       integration_result, błąd jest definiowany jako:
    %       integration_error(1,i) = abs(integration_result - reference_value),
    %       gdzie reference_value to wartość referencyjna całki.
    %
    %   Nt - wektor wierszowy zawierający liczby losowań, dla których obliczano
    %       wektor błędów całkowania integration_error.
    %
    %   ft_5 - gęstość funkcji prawdopodobieństwa dla n=5
    %
    %   [xr, yr] - tablice komórkowe zawierające informacje o wylosowanych punktach.
    %       Tablice te mają rozmiar [1, length(Nt)]. W komórkach xr{1,i} oraz yr{1,i}
    %       zawarte są współrzędne x oraz y wszystkich punktów zastosowanych
    %       do obliczenia całki przy losowaniu Nt(1,i) punktów.
    %
    %   yrmax - maksymalna dopuszczalna wartość współrzędnej y losowanych punktów

    reference_value = 0.0473612919396179; % wartość referencyjna całki

    Nt = 5:50:10^4;
    integration_error = zeros(1, length(Nt));

    mikron = 10;
    sigma = 3;

    f = @(t) (1/(sigma * sqrt(2*pi))) * exp(-((t - mikron)^2)/(2*sigma^2));
    ft_5 = f(5);
    yrmax = f(5.1);

    xr = cell(1, length(Nt));
    yr = cell(1, length(Nt));

    for i = 1:length(Nt)
        [integration, xr{1,i}, yr{1,i}] = monte_carlo_method_integral(f, 0, 5, Nt(i), yrmax);
        integration_error(i) = abs(integration - reference_value);
    end

    loglog(Nt,integration_error);
    xlabel('Liczba podprzedziałów całkowania');
    ylabel('Błąd całkowania');
    title('Błąd całkowania w zależności od liczby podprzedziałów całkowania');
    grid on;

    saveas(gcf, 'zadanie4.png');
end

function [result, xr, yr] = monte_carlo_method_integral(f, a, b, n, yrmax)
    xr = rand(1, n) * (b - a) + a;
    yr = rand(1, n) * yrmax;
    under_curve = yr <= arrayfun(f, xr);
    result = (b - a) * yrmax * mean(under_curve);
end 
