
function control = AD2D_Map_ActionIdx_to_Control( idx )
   if(idx == 1) 
        control = [0, 0];
   elseif(idx == 2)
        control = [1, 0];
   elseif(idx == 3)
        control = [0, 1];
   elseif(idx == 4)
        control = [1, 1];
   else
        error('Unknow control idx')
   end
end
