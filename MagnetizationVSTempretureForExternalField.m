% Parameters
numSpins = 100;            % Number of spins in the lattice
J = 1;                     % Interaction strength
numIterations = 1000;      % Number of iterations for spin updates
numTemperatureSteps = 100; % Number of temperature steps
fields = [0, 0.5, 1];      % External magnetic field values

% Temperature range
T_min = 0.1;
T_max = 5;
temperatures = linspace(T_min, T_max, numTemperatureSteps);

% Initialize matrix to store magnetization values for different fields
magnetization = zeros(numTemperatureSteps, length(fields));

% Loop over each magnetic field value
for f = 1:length(fields)
    h = fields(f); % Current magnetic field
    
    % Monte Carlo simulation for each temperature
    for t = 1:numTemperatureSteps
        T = temperatures(t);
        beta = 1 / T;
        
        % Initialize spins randomly
        spins = 2 * randi([0, 1], numSpins, 1) - 1;
        
        % Simulate over a number of iterations
        for iter = 1:numIterations
            % Randomly select a spin
            i = randi(numSpins);
            
            % Calculate the energy change if this spin is flipped
            deltaE = 2 * J * spins(i) * (spins(mod(i-2, numSpins) + 1) + spins(mod(i, numSpins) + 1)) + 2 * h * spins(i);
            
            % Metropolis criterion
            if deltaE < 0 || rand() < exp(-deltaE * beta)
                spins(i) = -spins(i); % Flip the spin
            end
        end
        
        % Calculate average magnetization at this temperature
        magnetization(t, f) = mean(spins);
    end
end

% Plotting the magnetization vs. temperature for different magnetic fields
figure;
hold on;
colors = {'b', 'r', 'g'}; % Colors for different magnetic field plots

for f = 1:length(fields)
    plot(temperatures, magnetization(:, f), 'Color', colors{f}, 'LineWidth', 1.5, 'DisplayName', ['h = ' num2str(fields(f))]);
end

xlabel('Temperature (T)');
ylabel('Magnetization (M)');
title('Magnetization vs. Temperature for Different Magnetic Fields');
legend('show');
grid on;
hold off;
