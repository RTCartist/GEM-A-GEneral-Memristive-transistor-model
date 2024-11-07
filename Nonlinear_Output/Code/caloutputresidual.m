% Given residuals
residuals = [0, 0.03884, 0.04205, 0.04461, 0.0454, 0.04412, 0.04028, 0.03367, 0.0228, 0.01207, ...
            -0.00664, -0.01381, -0.01569, -0.01871, -0.02369, -0.02473, -0.01951, -0.01683, ...
            -0.0151, -0.01128, -0.01246, -0.01125, -0.00655, -0.003, -0.0019, -0.000680836, ...
             0.00227, 0.00351, 0.00537, 0.00586, 0.00468, 0.00584, 0.00495];

% Mean Error
mean_error = mean(residuals);
fprintf('Mean Error: %.6f\n', mean_error);

% Mean Absolute Error (MAE)
mae = mean(abs(residuals));
fprintf('Mean Absolute Error (MAE): %.6f\n', mae);

% Root Mean Square Error (RMSE)
rmse = sqrt(mean(residuals.^2));
fprintf('Root Mean Square Error (RMSE): %.6f\n', rmse);

% Standard Deviation of Residuals
std_dev = std(residuals);
fprintf('Standard Deviation: %.6f\n', std_dev);
