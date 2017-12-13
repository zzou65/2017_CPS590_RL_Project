
function feature_vec = AD2D_Generate_Feature()
    AD2D_Global_Vars;
    if(If_Use_AddFeature == 1)
        feature_vec = [Current_S_Pos(:); Current_S_Mag];
    else
        feature_vec = [];
    end
    state = Image_Projection_Mat_for_Feature * Current_Sol(:);
    feature_vec = [feature_vec; state];
    

