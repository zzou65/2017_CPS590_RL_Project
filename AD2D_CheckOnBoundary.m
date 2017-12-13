
function if_on_bdry = AD2D_CheckOnBoundary( current_pos )
    if(abs(current_pos(1))<1e-10 || abs(current_pos(1)-1) < 1e-10)
        if_on_bdry = 1;
    elseif(abs(current_pos(2))<1e-10 || abs(current_pos(2)-1) < 1e-10)
        if_on_bdry = 1;
    else
        if_on_bdry = 0;
    end
end
