classdef GaussianFunctions
      % Contain some key functions based on the Gaussian assumption
      
      methods (Static)
            
            % Kalman Prediction
            function state_pred = predict(state, motion_model)
                  % Performs linear/nonlinear (Extended) Kalman prediction
                  % INPUT:
                  %     state.x: object state mean, (state dimension) * 1
                  %     state.P: object state covariance
                  % OUTPUT:
                  %     state_pred.x: predicted object state mean, (state dimension) * 1
                  %     state_pred.P: predicted object state covariance
                  
                  state_pred.x = motion_model.f(state.x);
                  state_pred.P = motion_model.F(state.x)*state.P*motion_model.F(state.x)'+motion_model.Q;
                  
            end
            
            
            
            % Kalman Update
            function state_upd = update(state_pred, z, measure_model)
                  % Performs linear/nonlinear (Extended) Kalman update step
                  % INPUT:
                  %     state_pred.x: predicted object state mean, (state dimension) * 1
                  %     state_pred.P: predicted object state covariance
                  %     z: measurement data
                  % OUTPUT:
                  %     state_upd.x: updated object state mean, (state dimension) * 1
                  %     state_upd.P: updated object state covariance
                  
                  H = measure_model.H(state_pred.x);  % measurement model Jacobian
                  S = H * state_pred.P * H' + measure_model.R;  % innovation covariance
                  S = (S + S') / 2;  % make sure matrix S is positive definite
                  K = (state_pred.P*H')/S;  % Kalman gain
                  
                  state_upd.x = state_pred.x + K * (z - measure_model.h(state_pred.x));
                  state_upd.P = (eye(length(state_pred.x)) - K * H) * state_pred.P;
                  
            end
            
            
            
            % Calculate the Predicted Likelihood
            function pl = predicted_likelihood(state_pred, z, measure_model)
                  % Calculates the predicted likelihood in logarithm domain which is used for weight calculation
                  % OUTPUT:
                  %     pl: predicted likelihood in logarithmic scale, (number of measurements) * 1
                  
                  if isempty(z)
                        pl = log(0);  % no measurement so zero weight
                  else
                        pl = zeros(length(z(1, :)), 1);
                        
                        Hmu = measure_model.h(state_pred.x);  % H * predicted mean
                        H = measure_model.H(state_pred.x);  % measurement model Jacobian
                        S = H * state_pred.P * H' + measure_model.R;  % innovation covariance
                        S = (S + S') / 2;  % make sure matrix S is positive definite
                        
                        for i = 1 : length(z(1, :))  % each measurement
                              pl(i) = log_mvnpdf(z(:, i), Hmu, S);  % the probability that given measurement belongs to proposed Gaussian
                        end
                  end
                  
            end
            
            
            
            % Approximate Gaussian mixture as a single Gaussian
            function state_merg = moment_match(w, states)
                  % Approximate a Gaussian mixture as a single Gaussian using moment matching
                  % INPUT:
                  %     w: normalized weights of Gaussian components in logarithm domain
                  %     states: (number of Gaussian components) * 1 struct array
                  %     states.x: Gaussian component mean, (state dimension) * 1
                  %     states.P: Gaussian component covariance
                  % OUTPUT:
                  %     state_merg.x: merged object state mean, (state dimension) * 1
                  %     state_merg.P: merged object state covariance
                  
                  if length(w) == 1
                        state_merg = states;
                  else
                        w = exp(w);  % convert back to decimal
                        state_merg.x = 0; state_merg.P = 0;
                        
                        for i = 1 : length(w)  % each Gaussian component
                              state_merg.x = state_merg.x + w(i)*states(i).x;
                        end
                        for i = 1 : length(w)
                              state_merg.P = state_merg.P + w(i) * states(i).P + w(i) * (state_merg.x - states(i).x) * (state_merg.x - states(i).x)';
                        end
                  end
                  
            end
            
      end
      
end