
function if_on_mid = AD2D_CheckOnMid( current_pos )
    if(abs(current_pos(1)-0.5) < 1e-10)
        if_on_mid = 1;
    else
        if_on_mid = 0;
    end
end 

