% left_or_right = 1/2
function V_Field = AD2D_V_Handle(coord, left_or_right, midcoord) 
    V_Field = zeros(size(coord));
    if left_or_right == 1
        idx = find( coord(:, 1) <= midcoord + 1e-10 );
    else
        idx = find( coord(:, 1) >= midcoord - 1e-10 );
    end
    V_Field(idx, 2) = 1;
end