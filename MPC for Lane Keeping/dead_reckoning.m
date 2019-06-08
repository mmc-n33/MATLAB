%% ================================================= %%
% This code was used to check if the developed model makes sense
% ==================================================== %


% ========= %
% Parameter %
% ========= %

% physical parameters
data = load('reference_data.mat');
rx = data.x; ry = data.y; % lane coordinates (m)
Tt = 80; % total time (s)
Ts = Tt / length(rx); % sampling time (s)
lr = 0.35; % rear tire to center distance (m)
% simulation parameters
n = 100;
% starting values for signals
xs = 6; ys = 7; psis = -pi/4;
vs = 1; betas = 0;


% ========= %
% Simulation %
% ========= %

% preallocations
x = [xs*ones(1, n+1) ; ys*ones(1, n+1) ; psis*ones(1, n+1)]; % [x-pos, y-pos, inertial angle]
u = [vs*ones(1, n+1) ; betas*ones(1, n+1)]; % [velocity, velocity angle]

for k = 1 : n
      
      if (2 <= k) && (k < 25)
            u(1, k:end) = u(1, k-1) + 0.02; % straight-line acceleration
      end
      if (25 <= k) && (k < 50)
            u(2, k:end) = u(2, k-1) + 0.01; % turn left
      end
      if (50 <= k) && (k < 75)
            u(2, k:end) = u(2, k-1) - 0.02; % reverse turn wilder
      end
      if 75 <= k
            u(1, k:end) = u(1, k-1) + 0.06; % fast & furious mode 
            u(2, k:end) = u(2, k-1) + 0.03; 
      end
      
      [A, B, Const] = model_linearization([x(3, k) ; u(:, k)], [Ts ; lr]);
      x(:, k+1) = A*x(:, k) + B*u(:, k) + Const;
      
end

% ==== %
% Plot %
% ==== %

set(0, 'defaultLineLineWidth', 1)
set(0, 'defaultLineMarkerSize', 8) % this follows a hierachy, that is "markersize" is a subclass of "line"
figure
plot(x(1, 1), x(2, 1), 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x(1, 25), x(2, 25), 'ws', 'MarkerFaceColor', [0 1 0]); hold on;
plot(x(1, 50), x(2, 50), 'ws', 'MarkerFaceColor', [0 0.7 0]); hold on;
plot(x(1, 75), x(2, 75), 'ws', 'MarkerFaceColor', [0 0.4 0]); hold on;
plot(x(1, end), x(2, end), 'bd', 'MarkerFaceColor', 'b'); hold on;
plot(x(1, :), x(2, :), 'k.-'); hold off;
xlabel('(m)'); ylabel('(m)'); title('xy'); legend('Starts Here', '25th', '50th', '75th', 'Ends Here');
legend('location', 'best');

figure
plot(1, x(3, 1), 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(25, x(3, 25), 'ws', 'MarkerFaceColor', [0 1 0]); hold on;
plot(50, x(3, 50), 'ws', 'MarkerFaceColor', [0 0.7 0]); hold on;
plot(75, x(3, 75), 'ws', 'MarkerFaceColor', [0 0.4 0]); hold on;
plot(length(x(3, :)), x(3, end), 'bd', 'MarkerFaceColor', 'b'); hold on;
plot(x(3, :), 'k.-'); hold off;
title('\psi'); legend('Starts Here', '25th', '50th', '75th', 'Ends Here');
legend('location', 'best');

figure
plot(1, u(2, 1), 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(25, u(2, 25), 'ws', 'MarkerFaceColor', [0 1 0]); hold on;
plot(50, u(2, 50), 'ws', 'MarkerFaceColor', [0 0.7 0]); hold on;
plot(75, u(2, 75), 'ws', 'MarkerFaceColor', [0 0.4 0]); hold on;
plot(length(u(2, :)), u(2, end), 'bd', 'MarkerFaceColor', 'b'); hold on;
plot(u(2, :), 'k.-'); hold off;
title('\beta'); legend('Starts Here', '25th', '50th', '75th', 'Ends Here');
legend('location', 'best');
