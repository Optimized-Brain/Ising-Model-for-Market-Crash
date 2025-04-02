% Parameters
numSpins = 100;            % Number of spins in the lattice
J = 1;                     % Interaction strength
numIterations = 1000;      % Number of iterations for spin updates
numFieldSteps = 100;       % Number of external field steps
temperatures = [1, 2.5, 4]; % Different temperature values

% External magnetic field range
h_min = -2;
h_max = 2;
fields = linspace(h_min, h_max, numFieldSteps);

% Initialize matrix to store magnetization values for different temperatures
magnetization = zeros(numFieldSteps, length(temperatures));

% Loop over each temperature
for tempIdx = 1:length(temperatures)
    T = temperatures(tempIdx); % Current temperature
    beta = 1 / T;
    
    % Monte Carlo simulation for each external field value
    for hIdx = 1:numFieldSteps
        h = fields(hIdx); % Current magnetic field
        
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
        
        % Calculate average magnetization at this external field
        magnetization(hIdx, tempIdx) = mean(spins);
    end
end

% Plotting the magnetization vs. external field for different temperatures
figure;
hold on;
colors = {'b', 'r', 'g'}; % Colors for different temperature plots

for tempIdx = 1:length(temperatures)
    plot(fields, magnetization(:, tempIdx), 'Color', colors{tempIdx}, 'LineWidth', 1.5, 'DisplayName', ['T = ' num2str(temperatures(tempIdx))]);
end

xlabel('External Magnetic Field (h)');
ylabel('Magnetization (M)');
title('Magnetization vs. External Field for Various Temperatures');
legend('show');
grid on;
hold off;
