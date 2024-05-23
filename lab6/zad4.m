load('energy.mat');

[country, source, degrees, y_original, y_yearly, y_approximation, mse, msek] = zadanie4(energy);

function [country, source, degrees, x_coarse, x_fine, y_original, y_yearly, y_approximation, mse, msek] = zadanie4(energy)
country = 'Poland';
source = 'Coal';
degrees = 1:4;
x_coarse = [];
x_fine = [];
y_original = [];
y_yearly = [];
y_approximation = [];

if isfield(energy, country) && isfield(energy.(country), source)
  dates = energy.(country).(source).Dates;
  y_original = energy.(country).(source).EnergyProduction;

  n_years = floor(length(y_original) / 12);
  y_cut = y_original(end-12*n_years+1:end);
  y4sum = reshape(y_cut, [12 n_years]);
  y_yearly = sum(y4sum, 1)';

  N = length(y_yearly);
  P = (N-1)*10 + 1;

  x_coarse = linspace(-1, 1, N)';
  x_fine = linspace(-1, 1, P)';

  y_approximation = cell(1, N-1);
  mse = zeros(N-1, 1);
  msek = zeros(N-2, 1);

  for i = 1:N-1
    p = my_polyfit(x_coarse, y_yearly, i);
    y_approx_coarse = polyval(p, x_coarse);
    y_approx_fine = polyval(p, x_fine);

    y_approximation{i} = y_approx_fine;

    mse(i) = mean((y_yearly - y_approx_coarse).^2);

    if i > 1
      msek(i-1) = mean((y_approximation{i} - y_approximation{i-1}).^2);
    end
  end

  degrees = ceil(linspace(1, N-1, 4));
  figure;

  subplot(3, 1, 1);
  plot(x_coarse, y_yearly, 'DisplayName', 'Dane roczne');
  hold on;
  for i = 1:length(degrees)
    plot(x_fine, y_approximation{degrees(i)}, 'DisplayName', ['Stopień ', num2str(degrees(i))]);
  end
  hold off;
  legend('show');
  title('Dane roczne i aproksymacje');
  xlabel('x');
  ylabel('y');

  subplot(3, 1, 2);
  semilogy(mse);
  title('Błąd średniokwadratowy (MSE) w funkcji stopnia wielomianu');
  xlabel('Stopień wielomianu');
  ylabel('MSE (skala logarytmiczna)');

  subplot(3, 1, 3);
  semilogy(msek);
  title('Błąd różnicowy między aproksymacjami różnych stopni');
  xlabel('Stopień wielomianu');
  ylabel('Błąd różnicowy (skala logarytmiczna)');

  saveas(gcf, 'zadanie4.png');

else
  disp(['Dane dla (country=', country, ') oraz (source=', source, ') nie są dostępne.']);
end
end

function p = my_polyfit(x, y, deg)
  X = zeros(length(x), deg+1);
  for i = 1:deg+1
    X(:,i) = x.^(deg+1-i);
  end
  p = (X'*X)\(X'*y);
end
