function main
    n_max = 200;
    a = 10;
    r_max = a/2; 
    
    [circles, index_number] = generate_circles(a, r_max, n_max);
    plot_circles(a, circles, index_number); 
end

function [circles, index_number] = generate_circles(a, r_max, n_max)
    % a - długość boku kwadratu
    % r_max - maksymalny promień okręgu
    % n_max - maksymalna liczba okręgów
    % circles - macierz zawierająca informacje o okręgach
    % index_number - liczba okręgów
    index_number = 193064;
    L1 = mod(index_number, 10);

    circles = zeros([n_max, 3]);
    for i = 1:n_max
        while true
            collides = 0;
            r = r_max*rand(1);
            x = a*rand(1);
            y = a*rand(1);

            for j = 1:i-1
                distance_x = x - circles(j, 1);
                distance_y = y - circles(j, 2);
                distance = sqrt(distance_x^2 + distance_y^2);
                if distance <= r + circles(j, 3)
                    collides = 1;
                end
            end
            if (x + r <= a) && (x - r >= 0) && (y + r <= a) && (y - r >= 0) && collides == 0
                break
            end
        end
        circles(i, 1) = x;
        circles(i, 2) = y;
        circles(i, 3) = r;
    end

    if mod(L1, 2) == 1
        circles = transpose(circles);
    end
end

function plot_circles(a, circles, index_number)
    % a - długość boku kwadratu
    % circles - macierz zawierająca informacje o okręgach
    % index_number - liczba okręgów
    figure
    hold on
    axis equal
    axis([0 a 0 a])
    L1 = mod(index_number, 10);
    if mod(L1, 2) == 0
        for i = 1:size(circles, 1)
            plot_circle(circles(i, 3), circles(i, 1), circles(i, 2));
        end
    else
        for i = 1:size(circles, 2)
            plot_circle(circles(3, i), circles(1, i), circles(2, i));
        end
    end

    hold off
    print -dpng 'zadanie1.png'
end