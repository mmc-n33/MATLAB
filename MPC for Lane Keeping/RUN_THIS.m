% ********************************************************** %
% This code designs an MPC controller and shows associated simulation %
% ......................................................................................................... %

% Major Problem: cannot find feasible solution to pass that sharp turn in a fast enough speed

%%
% ============ %
% PARAMETERS %
% ============ %

addpath(genpath('C:\Users\m8mao\Desktop\qqnbxl\Control\MPC\MPC Vehicle Simulation'));

% Physical parameters
data = load('reference_data.mat');
rx = data.x; ry = data.y; % lane coordinates (m)
rxo = data.xo; ryo = data.yo; % outer track coordinates (m)
rxi = data.xi; ryi = data.yi; % outer track coordinates (m)
Tt = 60; % total time (s)
Ts = 10*Tt / length(rx); % sampling time (s) (10* so that it is large enough to make more sense)
lr = 0.35; % rear tire to center distance (m)
xs = rx(1); ys = ry(1); psis = -pi; vs = 0.4; betas = 0; % specify vehicle's starting position

% MPC parameters
w = data.w; % track width (m)
vlim = 1; % velocity limit (m/s)
% vdlim = 0.2; % velocity difference limit (m/s)
betalim = atan(0.5*tan(40/180*pi)); % velocity angle limit (rad) (steering angle range is [-40 deg, 40 deg])
N = 7; % prediction horizon
Nc = 3; % control horizon
Q = eye(2); R = 0.1; Rd = 0.3*eye(2); P = -0.5; S = 100;


% ============= %
% OPTIMIZATION %
% ============= %

% Preallocation
% these are system states
x = [xs*ones(1, length(rx)+1) ; ... % x-pos
      ys*ones(1, length(rx)+1) ; ... % y-pos
      psis*ones(1, length(rx)+1)]; % vehicle angle
u = [vs*ones(1, length(rx)) ; ... % velocity
      betas*zeros(1, length(rx))]; % velocity angle (use beta instead of tire angle to avoid nonlinearity)
% these are prediction-horizon states
x_yal = sdpvar(3*ones(1, N+1), 1*ones(1, N+1)); % to distinguish x(k+1) and u(k), x_yal{1} will be skipped
u_yal = sdpvar(2*ones(1, N), 1*ones(1, N));
% s = sdpvar; % slack variable for tracking relaxation
% these are plot states
pred = ones(7, N, length(rx)-N+1); % store prediction results

options = sdpsettings('solver', 'gurobi', 'verbose', 0, 'debug', 0);
dibage = 0; % set to 1 to see "real-time" results for debugging purpose

tic;
for k = 1 : floor(length(rx)*1)
      
      % Propagate Yalmip constraints and objectives
      constraints = [];
      objective = 0;
      for j = 1 : N
            
            if k + j <= length(rx)
                  e_pos = [rx(k+j) - x_yal{j+1}(1) ; ry(k+j) - x_yal{j+1}(2)];
            else
                  e_pos = [rx(k+j-length(rx)) - x_yal{j+1}(1) ; ry(k+j-length(rx)) - x_yal{j+1}(2)];
            end
            
            objective = objective + e_pos' * Q * e_pos; % tracking
            %{
            objective = objective + P * u_yal{j}(1)^2; % wozaiganshijian
            objective = objective + S * s^2; % slack variable
            objective = obejctive +  P * (psi_s{j} - psi_s{j+1})^2; % heading angle stability
            %}
            
            if j <= Nc % stop updating u after control horizon
                  objective = objective + (u_yal{j+1} - u_yal{j})' * Rd * (u_yal{j+1} - u_yal{j}); % vehicle should be stable
                  %{
                  objective = objective + R * u_s{j}^2;  % control effort
                  %}
            else
                  u_yal{j} = u_yal{Nc};
            end
            
            if j == 1 % "manually" set model because don't have value(sdpvar)
                  [A, B, Const] = model_linearization([x(3, k) ; u(:, k)], [Ts ; lr]);
                  constraints = [constraints, x_yal{j+1} == A*x(:, k) + B*u_yal{j} + Const];
                  %{
                  constraints = [constraints, -vdlim <= u_yal{j}(1) - u(1, k), u_yal{j}(1) - u(1, k) <= vdlim];
                  %}
            else
                  [A, B, Const] = model_linearization([value(x_yal{j}(3)) ; value(u_yal{j-1})], [Ts ; lr]);
                  constraints = [constraints, x_yal{j+1} == A*x_yal{j} + B*u_yal{j} + Const];
                  %{
                  constraints = [constraints, -vdlim <= u_yal{j}(1) - u_yal{j-1}(1), u_yal{j}(1) - u_yal{j-1}(1) <= vdlim];
                  %}
            end
            
            if k + j <= length(rx)
                  constraints = [constraints, (rx(k+j) - x_yal{j+1}(1))^2 + (ry(k+j) - x_yal{j+1}(2))^2 <= (w / 2)^2];
            else
                  constraints = [constraints, ...
                        (rx(k+j-length(rx)) - x_yal{j+1}(1))^2 + (ry(k+j-length(rx)) - x_yal{j+1}(2))^2 <= (w / 2)^2];
            end
            constraints = [constraints, 0 <= u_yal{j}(1) <= vlim];
            constraints = [constraints, -betalim <= u_yal{j}(2) <= betalim];
            %{
            constraints = [constraints, 0 <= s <= 0.1];
            %}
            
            diagnostic = optimize(constraints, objective, options);
            
            if dibage == 1
                  if diagnostic.problem ~= 0
                        warning('Something went wrong !!!')
                        break
                  else
                        pred(:, i) = [value(x_yal{j+1}) ; value(u_yal{j}) ; rx(k+j) ; ry(k+j)];
                        i = i + 1;
                  end
                  
                  fprintf('tracking = %f\n', (rx(k+j) - value(x_yal{j+1}(1)))^2 + (ry(k+j) - value(x_yal{j+1}(2)))^2)
                  fprintf('x = %f\n', value(x_yal{j+1}(1)))
                  fprintf('y = %f\n', value(x_yal{j+1}(2)))
                  fprintf('psi = %f\n', value(x_yal{j+1}(3)))
                  fprintf('v = %f\n', value(u_yal{j}(1)))
                  fprintf('beta = %f\n', value(u_yal{j}(2)))
            else
                  if diagnostic.problem ~= 0
                        warning('Stray too far from the track !!!')
                  end
                  if k + j <= length(rx)
                        pred(:, j, k) = [value(x_yal{j+1}) ; value(u_yal{j}) ; rx(k+j) ; ry(k+j)];
                  else
                        pred(:, j, k) = [value(x_yal{j+1}) ; value(u_yal{j}) ; rx(k+j-length(rx)) ; ry(k+j-length(rx))];
                  end
            end
            
      end
      
      if dibage == 1
            % Prediction results plot
            figure
            plot(pred(1, :), pred(2, :), '.-'); title('xy'); hold on;
            plot(pred(6, :), pred(7, :), 'o-.'); hold off;
            figure
            plot(pred(3, :)); title('\psi');
            figure
            plot(pred(4, :)); title('v');
            figure
            plot(pred(5, :)); title('\beta');
      end
      
      % System update
      u(:, k) = value(u_yal{1});
      [A, B, Const] = model_linearization([x(3, k) ; u(:, k)], [Ts ; lr]);
      x(:, k+1) = A*x(:, k) + B*u(:, k) + Const;
      
      fprintf('Time elapsed: %0.8f\n', toc);
      
end


% ============= %
% MPC Result Plot %
% ============= %

figure('units', 'normalized', 'outerposition', [0 0 1 1]) % full screen
pos = [[0.05 0.15 0.5 0.8] ; [0.55, 0.5, 0.45, 0.45] ; [0.6, 0.15, 0.25, 0.3] ; [0.86, 0.22, 0.12, 0.15]];
plot_backgrounds(betalim, pos(1 : 2, :))
plot_animation(x, u, pred, [rx ; ry], [w, N, betalim], pos, 'gif')

rmpath(genpath('C:\Users\m8mao\Desktop\qqnbxl\Control\MPC\MPC Vehicle Simulation'));

% save('111.mat', 'x', 'u', 'pred')