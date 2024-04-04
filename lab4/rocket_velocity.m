function velocity_delta = rocket_velocity(t)
    % velocity_delta - różnica pomiędzy prędkością rakiety w czasie t oraz zadaną prędkością M
    % t - czas od rozpoczęcia lotu rakiety dla którego ma być wyznaczona prędkość rakiety
    M = 750; % [m/s]
    m0 = 150000;
    u = 2000;
    q = 2700;
    g = 1.622;
    
    velocity_delta = abs(u * log(m0 / (m0 - q * t)) - g * t) - M; 
end