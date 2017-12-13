
function AD2D_Initialize_PathToFollow( path )
    AD2D_Global_Vars;
    AD2D_ResetState();

    Current_S_Pos_Idx = path(1);
    Current_S_Pos = SPos(Current_S_Pos_Idx, :);

    If_Follow_Existing_Path = 1;
    Existing_Path = path;
    

