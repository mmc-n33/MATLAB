function measurement = generate_measurement(state, sensor_model, measure_model)
% Generate objects measurements and clutters measurements
% INPUT:
%       state: (number of total time step) * 1 cell array, (object state dimension) * (number of objects) each cell
% OUTPUT:
%       measurement: (number of total time step) * 1 cell array, (measurement dimension) * (number of measurements) each cell

measurement = cell(length(state), 1);

for i = 1 : length(state)  % each time
      % Object
      if ~isempty(state{i})  % if there is object
            good = (rand(length(state{i}(1, :)), 1) <= sensor_model.P_D)';  % > detection probability means object is not detected
            for j = 1 : length(good)
                  if good(j) == 1  % true
                        measurement{i} = [measurement{i} mvnrnd(measure_model.h(state{i}(:, j)), measure_model.R)'];
                  end
            end
      end
      
      % Clutter
      N_c = poissrnd(sensor_model.rate_c);  % number of clutters
      C = rand(measure_model.d, N_c) .* ...
            (sensor_model.range_c(:, 2) - sensor_model.range_c(:, 1)) + sensor_model.range_c(:, 1);  % positon in surveillance area
      
      % Overall measurements
      measurement{i} = [measurement{i} C]; 
end

end