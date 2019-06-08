classdef GenerateState
      % Generate object states data
      
      properties
            x_init  % initial state
            t_birth  % objects start to exist from this time
            t_death  % objects stop to exist after this time
            T  % total time / time step
      end
      
      methods (Static)
            
            % Initialization
            function self = init(x_init, t_birth, t_death, T)  % note that use "self" because Python, in MATLAB it can be any name
                  self.x_init = x_init;
                  self.t_birth = t_birth;  % (number of objects) * 1 array
                  self.t_death = t_death;
                  self.T = T;
            end
            
            
            
            % Generate Object States Data
            function x_data = generate_state(self, motion_model, if_noisy)
                  % INPUT:
                  %     motion_model: a specific motion model
                  %     if_noisy: 0 means no noise, 1 means yes noise
                  % OUTPUT:
                  %     x_data: (number of total time step) * 1 cell array, (object state dimension) * (number of objects) each cell
                  
                  x_data = cell(self.T, 1);
                  
                  if if_noisy == true
                        for i = 1 : length(self.t_birth)  % each object
                              x_pred = self.x_init(:, i);
                              x_data{1} = self.x_init(:, i);
                              for j = 2 : self.T  % each time step
                                    if j >= self.t_birth(i) && j <= min(self.t_death(i), self.T)  % "min" makes sure object does exist at this time
                                          % the first element in x_data is the prediction of what's after IC not IC itself
                                          x_data{j} = [x_data{j} mvnrnd(motion_model.f(x_pred), motion_model.Q)'];  % Gaussian
                                          x_pred = x_data{j}(:, end);
                                    end
                              end
                        end
                  else
                        for i = 1 : length(self.t_birth)  % each object
                              x_pred = self.x_init(:, i);
                              x_data{1} = self.x_init(:, i);
                              for j = 2 : self.T  % each time step
                                    if j >= self.t_birth(i) && j <= min(self.t_death(i), self.T)  % "min" makes sure object does exist at this time
                                          % the first element in x_data is the prediction of what's after IC not IC itself
                                          x_data{j} = [x_data{j} mvnrnd(motion_model.f(x_pred), zeros(length(self.x_init(:, 1))))'];  % Gaussian
                                          x_pred = x_data{j}(:, end);
                                    end
                              end
                        end
                  end

            end
            
      end
      
end
            
            
            
            