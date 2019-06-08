classdef SOT
      % Contains single-object-tracking algorithms 

      properties
            gating   % gating parameter
            reduction  % hypothesis reduction parameter
      end
      
      methods (Static)
            
            % Initialization
            function self = initialize(gating_size, w_min, M, merge_threshold) % note that use "self" because Python, in MATLAB it can be any name
                  % INPUT: 
                  %     w_min: allowed minimum hypothesis weight (prune)
                  %     M: allowed maximum number of hypotheses (cap)
                  %     merge_threshold: merging threshold (merge)

                  self.gating_size = gating_size;
                  self.reduction.w_min = log(w_min);
                  self.reduction.merge_threshold = merge_threshold;
                  self.reduction.M = M;
                  
            end
            
            
            
            % Nearest Neighbor
            function estimates = NN_filter(self, state, Z, motion_model, measure_model, sensor_model)
                  % Approximate using the nearest neighbor to current state
                  % INPUT: 
                  %     state.x: object state mean, (state dimension) * 1
                  %     state.P: object state covariance
                  %     Z: measurement data, 
                  %         (number of total time step) * 1 cell array, (measurement dimension) * (number of measurements) each cell
                  % OUTPUT:
                  %     estimates: estimated object states
                  %                         (number of total time step) * 1 cell array, (object state dimension) * (number of objects) each cell
                  
                  estimates = cell(length(Z), 1);
                  
                  for i = 1 : length(Z)  % each time
                        % Gating
                        z_in = HypothesesReduction.gating(state, Z{i}, measure_model, self.gating_size);

                        % Weight
                        pl = GaussianFunctions.predicted_likelihood(state, z_in, measure_model);
                        w = log(sensor_model.P_D) + pl - log(sensor_model.intensity_c);
                        [~, ind] = max(w);
                        
                        % Update
                        if w(ind) >= log(1 - sensor_model.P_D)  % compare with the "no detection" hypothesis
                              % use z_in(ind) which gives wrong result, stuck for a long time to fix it, add this comment for memory
                              state = GaussianFunctions.update(state, z_in(:, ind), measure_model);
                        end
                        estimates{i} = state.x;  % posterior
                        
                        % Predict
                        state = GaussianFunctions.predict(state, motion_model);
                  end
                  
            end
            
            
            
            % Probalistic Data Association
            function estimates = PDA_filter(self, state, Z, motion_model, measure_model, sensor_model)
                  % Approximate using probalistic data association
                  % INPUT: 
                  %     state.x: object state mean, (state dimension) * 1
                  %     state.P: object state covariance
                  %     Z: measurement data, 
                  %         (number of total time step) * 1 cell array, (measurement dimension) * (number of measurements) each cell
                  % OUTPUT:
                  %     estimates: estimated object states
                  %                         (number of total time step) * 1 cell array, (object state dimension) * (number of objects) each cell
                  
                  estimates = cell(length(Z), 1);
                  
                  for i = 1 : length(Z)  % each time
                        % Gating
                        z_in = HypothesesReduction.gating(state, Z{i}, measure_model, self.gating_size);
                        
                        % Weights
                        pl = GaussianFunctions.predicted_likelihood(state, z_in, measure_model);
                        w = [log(sensor_model.P_D) + pl - log(sensor_model.intensity_c) ; log(1 - sensor_model.P_D)];
                        w = normalize_log_weights(w);
                        
                        % Prune
                        [w_prune, ijk] = HypothesesReduction.prune(w, self.reduction.w_min);
                        w = normalize_log_weights(w_prune);
                        
                        % Update each Gaussian
                        state_upd = [];
                        if length(ijk) ~= 1  % = 1 means only the "no object detection" hypothesis
                              for j = 1 : length(ijk) - 1  % the last hypothesis is always "no object detection" which means no update
                                    state_upd = [state_upd ; GaussianFunctions.update(state, z_in(:, ijk(j)), measure_model)];
                              end
                        end
                        state = [state_upd ; state];
                        
                        % Merge Gaussians
                        state = GaussianFunctions.moment_match(w, state);
                        estimates{i} = state.x;  % posterior
                        
                        % Predict
                        state = GaussianFunctions.predict(state, motion_model);
                  end
                  
            end
            
            
            
            % Gaussian Sum
            function estimates = GS_filter(self, state, Z, motion_model, measure_model, sensor_model)
                  % Approximate using Gaussian sum filter
                  % INPUT: 
                  %     state.x: object state mean, (state dimension) * 1
                  %     state.P: object state covariance
                  %     Z: measurement data, 
                  %         (number of total time step) * 1 cell array, (measurement dimension) * (number of measurements) each cell
                  % OUTPUT:
                  %     estimates: estimated object states
                  %                         (number of total time step) * 1 cell array, (object state dimension) * (number of objects) each cell
                  
                  % missing?
                  
                  
            end
            
      end
      
end



