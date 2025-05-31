function W_el = WeightEl(matrix_el)
        
W_el = diag(sind(matrix_el(:,1))) .* diag(sind(matrix_el(:,1)));
end

