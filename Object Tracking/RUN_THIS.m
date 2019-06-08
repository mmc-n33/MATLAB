%% Define Models
% Create motion_model model
T = 1; sigmaV = 1; sigmaOmega = pi/180;
motion_model = MotionModel.ctmodel(T, sigmaV, sigmaOmega);

% Create measure_model model
sigma_r = 5; sigma_b = pi/180; s = [300 ; 400];
measure_model = MeasureModel.rangebearingmeasmodel(sigma_r, sigma_b, s);

% Create senor model
P_D = 0.7; rate_c = 60; range_c = [-1000 1000 ; -pi pi];
sensor_model = SensorModel.sensor_model(P_D, rate_c, range_c);

% %% Define Models
% % Create motion_model model
% T = 0.1; sigma_q = 5;
% motion_model = MotionModel.cvmodel(T, sigma_q);
% 
% % Create measure_model model
% sigma_r = 0.2;
% measure_model = MeasureModel.cvmeasmodel(sigma_r);
% 
% % Create senor model
% P_D = 0.95; rate_c = 10; range_c = [-50 50 ; -50 50];
% sensor_model = SensorModel.sensor_model(P_D, rate_c, range_c);
% 
% % Generate objects' states data
% T = 20; initial_state.x = [6 6 6 6]'; initial_state.P = eye(4);
% states = GenerateState.generate_state(GenerateState.init(initial_state.x, 1, T+1, T), motion_model, 1);
% 
% % Generate measure_model data
% measurements = generate_measurement(states, sensor_model, measure_model);

%% Generate Data
% Generate objects' states data
T = 100; initial_state.x = [0 0 10 0 pi/180]'; initial_state.P = diag([1 1 1 1*pi/180 1*pi/180].^2);
states = GenerateState.generate_state(GenerateState.init(initial_state.x, 1, T+1, T), motion_model, 1);

% Generate measure_model data
measurements = generate_measurement(states, sensor_model, measure_model);

%% Solve SOT and Plot
figure
for i = 1 : 4
      % SOT setting
      P_G = [0.8, 0.999, 0.999, 0.999];  % gating percentage
      gating_size = chi2inv(P_G, measure_model.d);  % chi-square distributed
      w_min = [1e-3, 1e-1, 1e-3, 1e-3]; M = [100, 100, 20, 100]; merging_threshold = [2 2 2 0.2];
      
      % SOT algorithms
      states_NN = SOT.NN_filter(SOT.initialize(gating_size(i), w_min(i), M(i), merging_threshold(i)), ...
            initial_state, measurements, motion_model, measure_model, sensor_model);
      states_PDA = SOT.PDA_filter(SOT.initialize(gating_size(i), w_min(i), M(i), merging_threshold(i)), ...
            initial_state, measurements, motion_model, measure_model, sensor_model);
      states_GS = SOT.GS_filter(SOT.initialize(gating_size(i), w_min(i), M(i), merging_threshold(i)), ...
            initial_state, measurements, motion_model, measure_model, sensor_model);
      
      % Extract estimates
      state = cell2mat(states');
      state_NN = cell2mat(states_NN'); state_PDA = cell2mat(states_PDA'); state_GS = cell2mat(states_GS');
      
      % Evaluate estimates
      error_NN = RMSE(states_NN, states); error_PDA = RMSE(states_PDA, states); error_GS = RMSE(states_GS, states);
      
      % Plot
      if i <=2 subplot(10, 2, [3, 5, 7, 9] + i - 1); else subplot(10, 2, [13, 15, 17, 19] + i - 3); end
      
      grid on; hold on;
      h2 = plot(state(1, :), state(2, :), 'g', 'linewidth', 4);
      h3 = plot(state_NN(1, :), state_NN(2, :), '-s', 'color', [0.960, 0.098, 0.921]);
      h4 = plot(state_PDA(1, :), state_PDA(2, :), '-d', 'color', [0.098, 0.682, 0.960]);
      h5 = plot(state_GS(1, :), state_GS(2, :), '-o', 'color', [0.960, 0.596, 0.098]);
      h1 = plot(initial_state.x(1), initial_state.x(2), 'ro', 'markersize', 7, 'markerfacecolor', 'r');
      if i == 4  % cannot use h1-5 for the overall legend because it overwrites the last legend, so need new lines
            ho2 = plot(state(1, :), state(2, :), 'g', 'linewidth', 2);
            ho3 = plot(state_NN(1, :), state_NN(2, :), '-s', 'color', [0.960, 0.098, 0.921]);
            ho4 = plot(state_PDA(1, :), state_PDA(2, :), '-d', 'color', [0.098, 0.682, 0.960]);
            ho5 = plot(state_GS(1, :), state_GS(2, :), '-o', 'color', [0.960, 0.596, 0.098]);
            ho1 = plot(initial_state.x(1), initial_state.x(2), 'ro', 'markersize', 7, 'markerfacecolor', 'r');
      end
      hold off;
      
      xlabel('Surveillance Area - x Coordinate (m)'); ylabel('Surveillance Area - y Coordinate (m)');
      title({['Gating % = ', num2str(P_G(i)), ', w_{min} = ', num2str(w_min(i)), ...
            ', H_{max} = ', num2str(M(i)), ', Merge Threshold = ', num2str(merging_threshold(i))]}, ...
            'fontname', 'Calibri');
      xl = xlim; yl = ylim; xr = 0.05*(xl(2) - xl(1)); yr = 0.05*(yl(2) - yl(1));
      xlim([xl(1) - xr, xl(2) + xr]); ylim([yl(1) - yr, yl(2) + yr]);
      
      set(gca, 'fontsize', 9);
      legend([h3, h4, h5], ['RMSE: ', num2str(error_NN, 4)], ['RMSE: ', num2str(error_PDA, 4)], ...
            ['RMSE: ', num2str(error_GS, 4)], 'Location', 'best', 'fontsize', 7)
      
end
lgd = subplot(10, 2, 1.5); axis(lgd, 'off');
pos = get(lgd, 'position');
t = title(['Motion Model: Constant Yaw Rate  |  Measure Model: Azimuth  |  Sensor Model: P^D = ', ...
      num2str(P_D), ', \lambda_c = ', num2str(rate_c)]);  % "overall" title
% t = title(['Motion Model: Constant Velocity  |  Measure Model: Constant Velocity  |  Sensor Model: P^D = ', ...
%       num2str(P_D), ', \lambda_c = ', num2str(rate_c)]);  % "overall" title
l = legend(lgd, [ho1, ho2, ho3, ho4, ho5], 'Start Here', 'Ground Truth', ...
      'Nearest Neighbor', 'Probalistic Data Association', 'Gaussian Sum', ...
      'fontsize', 7, 'orientation' ,'horizontal', 'box', 'off', 'position', pos);  % "overall" legend


