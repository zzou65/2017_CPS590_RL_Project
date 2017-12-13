
function img = AD2D_MakeImgFromState( state )
    AD2D_Global_Vars;
    temp = Image_Projection_Mat * state(:);
    temp = reshape(temp, [N_Pixels_Per_Dim, N_Pixels_Per_Dim]);
    temp = flipud(temp);
    
    %img = temp;
    
    temp = max( min(temp, state_hig), state_low );
    img = ( temp - state_low ) / ( state_hig - state_low ) * 1;
    
    %img = mat2gray(temp, [state_low, state_hig]);
end
    

