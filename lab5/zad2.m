[nodes_Chebyshev, V, V2, original_Runge, interpolated_Runge, interpolated_Runge_Chebyshev] = zadanie2();

function [nodes_Chebyshev, V, V2, original_Runge, interpolated_Runge, interpolated_Runge_Chebyshev] = zadanie2()
% nodes_Chebyshev - wektor wierszowy zawierający N=16 węzłów Czebyszewa drugiego rodzaju
% V - macierz Vandermonde obliczona dla 16 węzłów interpolacji rozmieszczonych równomiernie w przedziale [-1,1]
% V2 - macierz Vandermonde obliczona dla węzłów interpolacji zdefiniowanych w wektorze nodes_Chebyshev
% original_Runge - wektor wierszowy zawierający wartości funkcji Runge dla wektora x_fine=linspace(-1, 1, 1000)
% interpolated_Runge - wektor wierszowy wartości funkcji interpolującej określonej dla równomiernie rozmieszczonych węzłów interpolacji
% interpolated_Runge_Chebyshev - wektor wierszowy wartości funkcji interpolującej wyznaczonej
%       przy zastosowaniu 16 węzłów Czebyszewa zawartych w nodes_Chebyshev 
    N = 16;
    x_fine = linspace(-1, 1, 1000);
    nodes_Chebyshev = get_Chebyshev_nodes(N);

    V = vandermonde_matrix(N);
    V2 = vandermonde_matrix_chebyshev(N, nodes_Chebyshev);
    original_Runge = 1 ./ (1 + 25 * x_fine.^2);

    subplot(2,1,1);
    plot(x_fine, original_Runge, 'DisplayName', '1/(1+25x^2)', 'LineWidth', 2.0);
    hold on;
    x_coarse = linspace(-1,1,N);
    y_coarse = 1 ./ (1 + 25 * x_coarse.^2);
    c_runge = V \ y_coarse';
    interpolated_Runge = polyval(flipud(c_runge), x_fine);
    plot(x_coarse, y_coarse, 'o', 'DisplayName', 'węzły równomiernie rozmieszczone');
    plot(x_fine, interpolated_Runge, 'DisplayName', ['N = ', num2str(N)]);
    hold off

    xlabel('x');
    ylabel('f(x)');
    title('Interpolacja funkcji Rungego');
    legend('show');

    subplot(2,1,2);
    plot(x_fine, original_Runge, 'DisplayName', '1/(1+25x^2)', 'LineWidth', 2.0);
    hold on;
    y_coarse = 1 ./ (1 + 25 * nodes_Chebyshev.^2);
    c_runge = V2 \ y_coarse';
    interpolated_Runge_Chebyshev = polyval(flipud(c_runge), x_fine);
    plot(nodes_Chebyshev, y_coarse, 'o', 'DisplayName', 'węzły Czebyszewa');
    plot(x_fine, interpolated_Runge_Chebyshev, 'DisplayName', ['N = ', num2str(N)]);
    hold off

    xlabel('x');
    ylabel('f(x)');
    title('Interpolacja funkcji sin(2\pix)');
    legend('show');

    saveas(gcf, 'zadanie2.png');
end

function nodes = get_Chebyshev_nodes(N)
    % oblicza N węzłów Czebyszewa drugiego rodzaju
    nodes = cos((0:N-1) * pi / (N-1));
end

function V = vandermonde_matrix(N)
    % Generuje macierz Vandermonde dla N równomiernie rozmieszczonych w przedziale [-1, 1] węzłów interpolacji
    x_coarse = linspace(-1,1,N);
    V = zeros(N,N);
    for i = 1:N
        V(:,i) = x_coarse.^(i-1);
    end
end

function V = vandermonde_matrix_chebyshev(N, nodes)
    % Generuje macierz Vandermonde dla N węzłów interpolacji zdefiniowanych w wektorze nodes
    V = zeros(N,N);
    for i = 1:N
        V(:,i) = nodes.^(i-1);
    end
end
