function [ A ] = fix_matrix( A , denom_factor )
% Fix to ensure positive definiteness from Danaher et al
% Divide each off-diagonal element by sum of absolute values of
% off-diagonal elements in its row

p = size(A, 1);

for cur_row = 1:p
    cur_sum = sum(abs(A(cur_row, :))) - 1;
    if cur_sum ~= 0
      A(cur_row, :) = A(cur_row, :) / (denom_factor * cur_sum);
    end
    
    % Make sure diagonal entries are still 1
    A(cur_row, cur_row) = 1;
end

% Final matrix is average of matrix with its transpose
 A = (A + A') / 2;

end

