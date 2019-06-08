function err = RMSE(s1, s2)
% Calculates the root mean square error between 2 state sequences
%INPUT: 
%       s1, s2: (number of total time step) * 1 cell array, (state dimension) * 1 each cell

err = 0;
for i = 1 : length(s1)
      err = err + mean((s1{i} - s2{i}).^2);   % state dimension mean
end
err = sqrt(err / length(s1));  % sequence length mean

end