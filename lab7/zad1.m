[integration_error, Nt, ft_5, integral_1000] = zadanie1();

function [integration_error, Nt, ft_5, integral_1000] = zadanie1()
    % Numeryczne całkowanie metodą prostokątów.
    % Nt - wektor zawierający liczby podprzedziałów całkowania
    % integration_error - integration_error(1,i) zawiera błąd całkowania wyznaczony
    %   dla liczby podprzedziałów równej Nt(i). Zakładając, że obliczona wartość całki
    %   dla Nt(i) liczby podprzedziałów całkowania wyniosła integration_result,
    %   to integration_error(1,i) = abs(integration_result - reference_value),
    %   gdzie reference_value jest wartością referencyjną całki.
    % ft_5 - gęstość funkcji prawdopodobieństwa dla n=5
    % integral_1000 - całka od 0 do 5 funkcji gęstości prawdopodobieństwa
    %   dla 1000 podprzedziałów całkowania

    reference_value = 0.0473612919396179; % wartość referencyjna całki

    Nt = 5:50:10^4;
    integration_error = [];

    mikron = 10;
    sigma = 3;

    f = @(t) (1/(sigma * sqrt(2*pi))) * exp(-((t - mikron)^2)/(2*sigma^2));
    ft_5 = f(5);
    [integral_1000] = rectangle_method_integral(f, 0, 5, 1000);

    for i = 1:length(Nt)
        integration_result = rectangle_method_integral(f, 0, 5, Nt(i));
        integration_error(i) = abs(integration_result - reference_value);
    end

    loglog(Nt,integration_error);
    xlabel('Liczba podprzedziałów całkowania');
    ylabel('Błąd całkowania');
    title('Błąd całkowania w zależności od liczby podprzedziałów całkowania');
    grid on;

    saveas(gcf, 'zadanie1.png');
end

function [result] = rectangle_method_integral(f, a, b, n)
    h = (b - a) / n;
    x = a:h:b;
    result = zeros(1, n);

    for i = 1:n
        result(i) = f((x(i) + x(i+1)) / 2);
    end

    result = h * sum(result);
end
