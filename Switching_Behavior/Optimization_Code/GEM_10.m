% Main script for fitting procedure
%--------------------------------------------------------------------------
% This script performs switching behavior optimization using the GEM model
clear all;
close all;
clc;

% Define initial parameters to optimize
koff    =   0.1;
kon     =   -0.1;
boff    =  1;
bon     =  1;
input_parameters = [koff, kon, boff, bon];

% Define the lower and upper bounds for the parameters
lb  = [-10000, -10000, -10, -10]; % Lower bounds [k_off_lb, k_on_lb, b_off_lb, b_on_lb];
ub  = [10000, 10000, 10, 10];     % Upper bounds [k_off_ub, k_on_ub, b_off_ub, b_on_ub];

% Read reference data from Excel file
data = readmatrix('..\Simulated_Data\Reference10\switching.xlsx');  % Replace with file name

v_ref = data(:, 1); % Column for voltage
v_ds  = 0.10;       % Drain-source bias voltage
i_ref = data(:, 2); % Column for current
t_ref = data(:, 3); % Column for time

% Perform optimization using two algorithms
optres1 = simulannealbnd(@(input_parameters) mse(input_parameters, v_ref, i_ref, t_ref), input_parameters, lb, ub);
optres2 = fminunc(@(input_parameters) mse(input_parameters, v_ref, i_ref, t_ref), input_parameters);

% Visual interpretation of optimization results
D = 1; % Device length
x_init = 0; % Initial value for state variable [0,D]
iv = 1; % I-V relationship, linear = 0, exponential = 1

% Run the GEM model with optimized parameters
[v_model1, i_model1, ~] = GEM_model(v_ref, v_ds, t_ref, x_init, iv, optres1);
[v_model2, i_model2, ~] = GEM_model(v_ref, v_ds, t_ref, x_init, iv, optres2);

% Calculate and display the MSE for each optimization
disp('MSE for Simulated Annealing:'), disp(mse(optres1, v_ref, i_ref, t_ref));
disp('MSE for Gradient Descent:'), disp(mse(optres2, v_ref, i_ref, t_ref));

% Plot the results
figure;
subplot(2,1,1);
hold on;
plot(t_ref, i_ref, 'r');
plot(t_ref, i_model1, 'b');
title('Simulated Annealing');

subplot(2,1,2);
hold on;
plot(t_ref, i_ref, 'r');
plot(t_ref, i_model2, 'b');
title('Gradient Descent');
hold off;

%--------------------------------------------------------------------------
% MSE Function for Optimization
%--------------------------------------------------------------------------
function out = mse(input_parameters, v_ref, i_ref, t_ref)
    % Model parameters
    x_init = 0; % Initial value for state variable [0,D]
    iv = 1; % I-V relationship, linear = 0, exponential = 1
    v_ds = 0.1; % Drain-source bias voltage

    % Run the GEM model with the given parameters
    [~, i_model, ~] = GEM_model(v_ref, v_ds, t_ref, x_init, iv, input_parameters);

    % Calculate RMSE as the fitting criterion
    N = length(t_ref);
    out = sqrt(sum((i_model - i_ref).^2) / sum(i_ref.^2) / N);
end

%--------------------------------------------------------------------------
% GEM Model Function
%--------------------------------------------------------------------------
function [v, i, x] = GEM_model(v_input, v_ds, t_vec, x_init, iv, input_parameters)
    % Fill model parameters 
    Roff = 16;
    Ron = 0.024;
    voff = -1.9; % Must be negative
    von = 8; % Must be positive
    koff = input_parameters(1);
    kon = input_parameters(2);
    aoff = 1;
    aon = 1;
    soff = 1;
    son = 1;
    boff = input_parameters(3);
    bon = input_parameters(4);
    D = 1;
    xoff = 0;
    xon = 1; % 1 corresponds to low-resistance state
    
    % Initialize variables
    i = zeros(length(t_vec), 1);
    v = v_input;
    lambda = reallog(Ron / Roff);

    % Simulation loop
    for j = 1:length(t_vec)
        if j == 1
            x(j) = x_init;
            if iv == 0
                i(j) = v_ds / (Roff + ((Ron - Roff) / (xon - xoff)) * (x(j) - xoff));
            elseif iv == 1
                i(j) = v_ds / (Roff * exp(lambda * (x(j) - xoff) / (xon - xoff)));
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
                i(j) = v_ds / (Roff * exp(lambda * (x(j) - xoff) / (xon - xoff)));
            end
        end
    end
end
