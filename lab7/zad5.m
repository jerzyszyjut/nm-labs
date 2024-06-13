[lake_volume, x, y, z, zmin] = zadanie5();

function [lake_volume, x, y, z, zmin] = zadanie5()
    % Funkcja zadanie5 wyznacza objętość jeziora metodą Monte Carlo.
    %
    %   lake_volume - objętość jeziora wyznaczona metodą Monte Carlo
    %
    %   x - wektor wierszowy, który zawiera współrzędne x wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   y - wektor wierszowy, który zawiera współrzędne y wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   z - wektor wierszowy, który zawiera współrzędne z wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   zmin - minimalna dopuszczalna wartość współrzędnej z losowanych punktów
    N = 1e6;
    zmin = -45;
    
    x_max = 100;
    y_max = 100;
    z_max = 0;

    x = x_max*rand(1,N); % [m]
    y = y_max*rand(1,N); % [m]
    z = zmin + (z_max - zmin)*rand(1,N);

    lake_depth = get_lake_depth(x, y);
    
    points_in_the_lake = z > lake_depth;
    lake_volume = x_max*y_max*(z_max - zmin);
    points_inside_count = sum(points_in_the_lake);
    lake_volume = lake_volume * points_inside_count/N;
end

function f = get_lake_depth( x, y )
    % Zwraca aproksymację głębokości jeziora w punkcie (x,y).
    % Zmienna f określa głębokość jeziora w metrach.
    x = x/100;
    y = y/100;

    w1=5; a1=30; x1=0.2; y1=0.3;
    w2=6; a2=40; x2=0.6; y2=1.0;
    w3=7; a3=50; x3=0.9; y3=0.5;
    w4=6; a4=40; x4=0.6; y4=0.1;
    w5=7; a5=70; x5=0.1; y5=0.95;
    w6=-5; a6=10; x6=0.5; y6=0.5;
    f=10 + w1*exp(-a1*((x-x1).^2+(y-y1).^2))+w2*exp(-a2*((x-x2).^2+(y-y2).^2))+w3*exp(-a3*((x-x3).^2+(y-y3).^2))+w4*exp(-a4*((x-x4).^2+(y-y4).^2))+w5*exp(-a5*((x-x5).^2+(y-y5).^2))+w6*exp(-a6*((x-x6).^2+(y-y6).^2));
    f = f - 9.5;
    f = f * 10;
    f(f>0) = 0;
end
