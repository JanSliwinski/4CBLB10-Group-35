function PlotThermoProperties(T, U, S, H, species_names)
    figure;
    subplot(3,1,1);
    plot(T, U);
    xlabel('Temperature (K)');
    ylabel('Internal Energy (J/mol)');
    legend(species_names);
    title('Internal Energy vs. Temperature');
    grid on;
    
    subplot(3,1,2);
    plot(T, H);
    xlabel('Temperature (K)');
    ylabel('Enthalpy (J/mol)');
    legend(species_names);
    title('Enthalpy vs. Temperature');
    grid on;
    
    subplot(3,1,3);
    plot(T, S);
    xlabel('Temperature (K)');
    ylabel('Entropy (J/mol·K)');
    legend(species_names);
    title('Entropy vs. Temperature');
    grid on;
end
