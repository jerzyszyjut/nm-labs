load('energy.mat');

[country, source, degrees, y_original, y_movmean, y_approximation, mse] = zadanie2(energy);

function [country, source, degrees, y_original, y_movmean, y_approximation, mse] = zadanie2(energy)
% Głównym celem tej funkcji jest wyznaczenie aproksymacji wygładzonych danych o produkcji energii elektrycznej w wybranym kraju i z wybranego źródła energii.
% Wygładzenie danych wykonane jest poprzez wywołanie funkcji movmean.
% Wybór kraju i źródła energii należy określić poprzez nadanie w tej funkcji wartości zmiennym typu string: country, source.
% Dopuszczalne wartości tych zmiennych można sprawdzić poprzez sprawdzenie zawartości struktury energy zapisanej w pliku energy.mat.
% 
% energy - struktura danych wczytana z pliku energy.mat
% country - [String] nazwa kraju
% source  - [String] źródło energii
% degrees - wektor zawierający cztery stopnie wielomianu dla których wyznaczono aproksymację
% y_original - dane wejściowe, czyli pomiary produkcji energii zawarte w wektorze energy.(country).(source).EnergyProduction
% y_movmean - średnia 12-miesięczna danych wejściowych, y_movmean = movmean(y_original,[11,0]);
% y_approximation - tablica komórkowa przechowująca cztery wartości funkcji aproksymującej wygładzone dane wejściowe. y_approximation stanowi aproksymację stopnia degrees(i).
% mse - wektor o rozmiarze 4x1: mse(i) zawiera wartość błędu średniokwadratowego obliczonego dla aproksymacji stopnia degrees(i).

country = 'Poland';
source = 'Coal';
degrees = [10, 20, 30, 40];
y_original = [];
y_approximation = [];
mse = [];

% Sprawdzenie dostępności danych
if isfield(energy, country) && isfield(energy.(country), source)
  % Przygotowanie danych do aproksymacji
  dates = energy.(country).(source).Dates;
  y_original = energy.(country).(source).EnergyProduction;
  y_movmean = movmean(y_original,[11,0]);

  x = linspace(-1,1,length(y_original))';

  % Pętla po wielomianach różnych stopni
  for i = 1:length(degrees)
    p = polyfit(x, y_movmean, degrees(i));
    y_approximation{i} = polyval(p, x);
    mse(i) = mean((y_movmean - y_approximation{i}).^2);
  end

  subplot(2,1,1);
  plot(dates, y_original, 'b', 'DisplayName', 'Dane oryginalne');
  hold on;
  plot(dates, y_movmean, 'r', 'DisplayName', 'Dane wygładzone');
  for i = 1:length(degrees)
    plot(dates, y_approximation{i}, 'DisplayName', ['Stopień wielomianu: ', num2str(degrees(i))]);
  end
  lgd = legend;
  lgd.NumColumns = 2;
  lgd.Location = 'southwest';
  lgd.FontSize = 6;
  legend('show');

  title('Aproksymacje i dane oryginalne');
  xlabel('Data');
  ylabel('Produkcja energii [GWh]');
  hold off;

  subplot(2,1,2);
  bar(mse);
  set(gca, 'XTickLabel', degrees);
  title('Błąd średniokwadratowy');
  xlabel('Stopień wielomianu');
  ylabel('MSE');
  grid on;

  saveas(gcf, 'zadanie2.png');

else
  disp(['Dane dla (country=', country, ') oraz (source=', source, ') nie są dostępne.']);
end

end
