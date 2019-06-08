classdef SensorModel
      % Creates a sensor model object
      
      methods (Static)
            
            % Calculate Sensor Model Parameters
            function obj = sensor_model(P_D, rate_c, range_c)
                  % INPUT:
                  %     P_D: probability that object is detected
                  %     rate_c: (lambda_c)^bar
                  %     range_c: range of clutter surveillance area, 1D: [xmin xmax], 2D: [xmin xmax ; ymin ymax]
                  
                  if length(range_c(:, 1)) > 1
                        area = (range_c(1, 2) - range_c(1, 1)) * (range_c(2, 2) - range_c(2, 1));
                  else
                        area = range_c(2) - range_c(1);
                  end

                  obj.intensity_c = rate_c / area;  % clutter intensity
                  obj.P_D = P_D;
                  obj.rate_c = rate_c;
                  obj.range_c = range_c;  
            end
            
      end
      
end