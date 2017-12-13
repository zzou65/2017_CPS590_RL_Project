
function AD2D_ResetState()
    AD2D_Global_Vars;

    Current_Time = 0;

    Current_Sol_S = zeros(NNode, 1);
    Current_Sol_C = zeros(NNode, 1);
    Current_Sol = Current_Sol_S + Current_Sol_C;

    Current_S_Pos = [1/4, 1/2];
    Current_S_Pos_Idx = AD2D_AsRdWalk([], SPos, stepSize);
    Previous_S_Pos = [-1, -1];
    Previous_S_Pos_Idx = -1;

    Current_S_Mag = 1.0;

    If_Follow_Existing_Path = 0;
    Existing_Path = [];


