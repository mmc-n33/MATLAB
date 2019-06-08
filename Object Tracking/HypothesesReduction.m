classdef HypothesesReduction
      % Contains methods to reduce number of hypotheses to simplify computation
      
      methods (Static)
            
            % Gating
            function z_in = gating(state_pred, z, measure_model, gating_size)
                  % Performs an ellipsoidal gating to eliminate measurements outside the gate
                  % INPUT:
                  %     state_pred.x: predicted object state mean, (state dimension) * 1
                  %     state_pred.P: predicted object state covariance
                  %     z: measurement data
                  %OUTPUT:
                  %     z_in: measurements in the gate
                  
                  Hmu = measure_model.h(state_pred.x);  % H * predicted mean
                  H = measure_model.H(state_pred.x);  % measurement model Jacobian
                  S = H * state_pred.P * H' + measure_model.R;  % innovation covariance
                  S = (S + S') / 2;  % make sure matrix S is positive definite
                  
                  z_in = double.empty(length(z(:, 1)), 0);
                  if ~isempty(z)
                        for i = 1 : length(z(1, :))  % each measurement
                              if (z(:, i) - Hmu)' / S * (z(:, i) - Hmu) <= gating_size
                                    z_in = [z_in z(:, i)];
                              end
                        end
                  end

            end
            
            
            
            % Prune
            function [w_new, ind] = prune(w, threshold)
                  % Prunes hypotheses with small weights
                  % OUTPUT:
                  % ind: the indices of hypotheses that are kept
                  
                  del = []; ind = [];
                  
                  for i = 1 : length(w)  % each hypothesis
                        if w(i) < threshold
                              del = [del i];  % cannot delete right here because it affects the indices of rest elements
                        else
                              ind = [ind i];
                        end
                  end
                  
                  w_new = w;
                  w_new(del) = [];
                  
            end
            
            
            
            % Cap
            function w_new = cap(w, M)
                  % Keeps M hypotheses with the highest weights and discard the rest
                  
                  [~, ind] = sort(w, 'descend');
                  ind_del = ind(M+1 : end);
                  
                  w_new = w;
                  w_new(ind_del) = [];

            end
            
            
            
            % Merge
            function [w_hat, states_hat] = mixture_reduction(w, states, threshold)
                  % Merges hypotheses within small Mahalanobis distance using a greedy method
                  % INPUT:
                  %     w: normalized weights of Gaussian components in logarithm domain
                  %     states: (number of Gaussian components) * 1 struct array
                  %     states.x: Gaussian component mean, (state dimension) * 1
                  %     states.P: Gaussian component covariance
                  %     threshold: discard if exceeding this threshold
                  % OUTPUT:
                  %     w_hat: normalized weights of remained Gaussian components in logarithmic scale
                  %     states_hat.x: remained Gaussian components means,
                  %                              (state dimension) * (number of remained Gaussian components)
                  %     states_hat.P: remained Gaussian components covariances,
                  %                              (state dimension) * (state dimension) * (number of remained Gaussian components)
                  
                  if length(w) == 1
                        w_hat = w;
                        states_hat = states;
                  else
                        for i = 1 : length(states)
                              all = [1 : length(states)];  % index set for all Gaussian components
                              
                              while ~isempty(all)
                                    good = [];
                                    [~, j] = max(w);
                                    
                                    for k = 1 : length(states)
                                          % Find all close components based on Mahalanobis distance
                                          if diag((states(j).x - states(k).x)' .* (states(j).P \ (states(j).x - states(k).x))) < threshold
                                                good = [good k];
                                          end
                                    end
                                    
                                    % Merge components
                                    [w_log, w_hat(i)] = normalize_log_weights(w(good));
                                    temp = GaussianFunctions.moment_match(w_log, states(good));
                                    states_hat.x = [states_hat.x temp]; states_hat.P = [states_hat.P temp];
                                    
                                    % Remove merged components from index set
                                    all(all == good) = [];
                              end
                        end
                  end
                  
            end
            
            
      end
      
end