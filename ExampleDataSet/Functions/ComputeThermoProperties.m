function [U, S, H] = ComputeThermoProperties(T, species_indices, SpS)
    nT = length(T);
    nSpecies = length(species_indices);
    
    U = zeros(nT, nSpecies);
    S = zeros(nT, nSpecies);
    H = zeros(nT, nSpecies);
    
    for i = 1:nSpecies
        species = SpS(species_indices(i));
        [u, s, h] = SpeciesThermo(T, species);
        U(:, i) = u;
        S(:, i) = s;
        H(:, i) = h;
    end
end

function [u, s, h] = SpeciesThermo(T, species)
    R = 8.314;  % J/mol·K
    nT = length(T);
    u = zeros(nT, 1);
    s = zeros(nT, 1);
    h = zeros(nT, 1);

    % Assuming only one temperature range
    T_low = species.Ts;  % Starting temperature
    T_high = species.Ts; % Ending temperature (since only one value is provided)
    a = species.Pol(1, :);  % Polynomial coefficients

    % Check that the coefficients are of length 7 (NASA polynomials)
    if length(a) < 7
        error('Polynomial coefficients must be of length 7.');
    end

    for i = 1:nT
        Ti = T(i);
        t = Ti / 1000;  % Scale temperature

        % NASA polynomials calculations
        cp_R = a(1) + a(2)*t + a(3)*t^2 + a(4)*t^3 + a(5)*t^4;
        h_RT = a(1) + a(2)*t/2 + a(3)*t^2/3 + a(4)*t^3/4 + a(5)*t^4/5 + a(6)/t;
        s_R = a(1)*log(t) + a(2)*t + a(3)*t^2/2 + a(4)*t^3/3 + a(5)*t^4/4 + a(7);

        h(i) = h_RT * R * Ti;
        s(i) = s_R * R;
        u(i) = h(i) - R * Ti;
    end
end
