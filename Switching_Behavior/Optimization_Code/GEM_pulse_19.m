% Main script for fitting procedure
%--------------------------------------------------------------------------
% This script performs switching behavior optimization using the GEM model
clear all;
close all;
clc;

% Define initial parameters to optimize
koff    =   -657.9194;
kon     =   2323.6;
boff    =  17.27;
bon     =  12.25;
input_parameters = [koff, kon, boff, bon];

% Define the lower and upper bounds for the parameters
lb  = [-1000, -1000, -200, -200]; % Lower bounds [k_off_lb, k_on_lb, b_off_lb, b_on_lb];
ub  = [1000, 1000, 200, 200];     % Upper bounds [k_off_ub, k_on_ub, b_off_ub, b_on_ub];

% Define pulse sequence
% Example: 50 pulses of 0.2V with 1ms duration, followed by 20 pulses of -0.1V with 2ms duration
data = readmatrix('..\Simulated_Data\Reference19\voltage.csv');
i_ref = data(:, 3);
% pulse_sequence = [2 * ones(20, 1); -1.4 * ones(20, 1)];
pulse_sequence = data(:, 2);
pulse_duration = [500e-6 * ones(8, 1); 500e-6 * ones(4, 1)];
t_ref = [0; cumsum(pulse_duration)]; % Adding an initial 0 time point
v_ds  = 0.1; 

% i_ref = randn(length(pulse_sequence), 1); % Random current values for testing
% Visual interpretation of optimization results
D = 1; % Device length
% x_init = 0.379; % Initial value for state variable [0,D] 
x_init = 0; % Initial value for state variable [0,D]
iv = 1; % I-V relationship, linear = 0, exponential = 1

% Perform optimization using two algorithms

optres1 = simulannealbnd(@(input_parameters) mse(input_parameters, pulse_sequence, i_ref, t_ref, x_init, iv, v_ds), input_parameters, lb, ub);
optres2 = fminunc(@(optres1) mse(optres1, pulse_sequence, i_ref, t_ref, x_init, iv, v_ds), optres1);

% Run the GEM model with optimized parameters
[v_model1, i_model1, ~] = GEM_model(pulse_sequence, v_ds, t_ref, x_init, iv, optres1);
[v_model2, i_model2, ~] = GEM_model(pulse_sequence, v_ds, t_ref, x_init, iv, optres2);

% Calculate and display the MSE for each optimization
disp('MSE for Simulated Annealing:'), disp(mse(optres1, pulse_sequence, i_ref, t_ref, x_init, iv, v_ds));
disp('MSE for Gradient Descent:'), disp(mse(optres2, pulse_sequence, i_ref, t_ref, x_init, iv, v_ds));

% Plot the results
figure;
subplot(2,1,1);
hold on;
plot(t_ref(2:end), i_ref, 'r'); % Exclude the initial 0 time point from plot
plot(t_ref(2:end), i_model1, 'b');
title('Simulated Annealing');

subplot(2,1,2);
hold on;
plot(t_ref(2:end), i_ref, 'r');
plot(t_ref(2:end), i_model2, 'b');
title('Gradient Descent');
hold off;

%--------------------------------------------------------------------------
% MSE Function for Optimization
%--------------------------------------------------------------------------
function out = mse(input_parameters, v_pulses, i_ref, t_ref, x_init, iv, v_ds)
    % Model parameters
    % x_init = 0; % Initial value for state variable [0,D]
    % iv = 2; % I-V relationship, linear = 0, exponential = 1
    % v_ds = 1; % Drain-source bias voltage

    % Run the GEM model with the given parameters
    [~, i_model, ~] = GEM_model(v_pulses, v_ds, t_ref, x_init, iv, input_parameters);

    % Calculate RMSE as the fitting criterion
    N = length(t_ref) - 1; % Exclude initial 0 time point
    out = sqrt(sum((i_model - i_ref).^2) / sum(i_ref.^2) / N);
end

%--------------------------------------------------------------------------
% GEM Model Function
%--------------------------------------------------------------------------
function [v, i, x] = GEM_model(v_pulses, v_ds, t_vec, x_init, iv, input_parameters)
    % Fill model parameters 
    % Roff = 36363636.36;
    % Ron = 57012.54276;
    Roff = 6.8e6;
    Ron  = 1.6e4;
    voff = -0.1; % Must be negative
    von = 2; % Must be positive
    koff = input_parameters(1);
    kon = input_parameters(2);
    aoff = 1;
    aon = 1;
    soff = 0.2;
    son = 0.2;
    boff = input_parameters(3);
    bon = input_parameters(4);
    D = 1;
    xoff = 0;
    xon = 1; % 1 corresponds to low-resistance state
    
    % Initialize variables
    i = zeros(length(t_vec)-1, 1); % Length matches the pulse sequence
    v = v_pulses;
    lambda = reallog(Roff / Ron);
    % beta = exp(Ron/Roff);
    % Simulation loop
    for j = 1:length(t_vec)-1
        if j == 1
            x(j) = x_init;
            if iv == 0
                i(j) = v_ds / (Roff + ((Ron - Roff) / (xon - xoff)) * (x(j) - xoff));
            elseif iv == 1
                i(j) = v_ds / (Ron / exp(-lambda * (x(j) - xon) / (xoff - xon)));
            elseif iv == 2
                i(j) = v_ds / (Ron * Roff / ((Roff - Ron) * reallog((exp(1) - 1) * (x(j) - xoff) / (xon - xoff) + 1) + Ron));
            end
        else
            dt = t_vec(j) - t_vec(j - 1);
            if v(j) <= voff
                dxdt(j) = (koff * ((v(j) / voff) - 1)^aoff) * (1 - x(j-1) / (xon - xoff) * soff)^boff;
                x(j) = x(j-1) + dxdt(j) * dt;
            elseif v(j) >= von
                dxdt(j) = (kon * ((v(j) / von) - 1)^aon) * (1 - x(j-1) / (xon - xoff) * son)^bon;
                x(j) = x(j-1) + dxdt(j) * dt;
            else
                dxdt(j) = 0;
                x(j) = x(j-1);
            end

            if x(j) < 0
                x(j) = 0;
                dxdt(j) = 0;
            elseif x(j) > D
                x(j) = D;
                dxdt(j) = 0;
            end
    
            if iv == 0
                i(j) = v_ds / (Roff + ((Ron - Roff) / (xon - xoff)) * (x(j) - xoff));
            elseif iv == 1
                i(j) = v_ds / (Ron / exp(-lambda * (x(j) - xon) / (xoff - xon)));
            elseif iv == 2
                i(j) = v_ds / (Ron * Roff / ((Roff - Ron) * reallog((exp(1) - 1) * (x(j) - xoff) / (xon - xoff) + 1) + Ron));
            end
        end
    end
end
