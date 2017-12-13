
function idx = AD2D_Map_Control_to_ActionIdx( control )
   if(control(1) == 0 && control(2) == 0) 
        idx = 1;
   elseif(control(1) == 1 && control(2) == 0)
        idx = 2;
   elseif(control(1) == 0 && control(2) == 1)
        idx = 3;
   elseif(control(1) == 1 && control(2) == 1)
        idx = 4;
   else
        error('Unknow control')
   end
end
