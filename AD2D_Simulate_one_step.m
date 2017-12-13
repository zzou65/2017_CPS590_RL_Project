% control [0, 0], [1, 0], [0, 1], [1, 1]

function [reward, scaled_reward] = AD2D_Simulate_one_step( current_control )
    AD2D_Global_Vars;
    
    if_on_mid = AD2D_CheckOnMid( Current_S_Pos );
    if (if_on_mid == 1 && current_control(1) == 1 && current_control(2) == 1)
        Current_S_Mag = 0.5 * Current_S_Mag;
    end
    if_on_bdry = AD2D_CheckOnBoundary( Current_S_Pos );
    if (if_on_bdry == 1 && Current_S_Pos(1) < 0.5 - 1e-10 && current_control(1) == 1 )
        Current_S_Mag = 0.5 * Current_S_Mag;
    end
    if (if_on_bdry == 1 && Current_S_Pos(1) > 0.5 + 1e-10 && current_control(2) == 1 )
        Current_S_Mag = 0.5 * Current_S_Mag;
    end
    
    [MatS, MatC, FSVec0, FCVec0] = Obtain_System( current_control );
    MatS = -MatS;
    MatC = -MatC;

    Previous_S_Pos_Idx = Current_S_Pos_Idx;
    Previous_S_Pos = Current_S_Pos;

    Current_Time = Current_Time + 1;

    if(If_Follow_Existing_Path ~= 1)
        Current_S_Pos_Idx = AD2D_AsRdWalk(Current_S_Pos_Idx, SPos, stepSize);
    else
        Current_S_Pos_Idx = Existing_Path(round(Current_Time + 1));
    end

    if Current_S_Pos_Idx <= 0
        Current_S_Pos = [-1, -1];
    else
        Current_S_Pos = SPos(Current_S_Pos_Idx, :);
    end

%%%%%%first time hit boundary decrease mag
%     if_on_bdry_0 = AD2D_CheckOnBoundary( Previous_S_Pos );
%     if_on_bdry_1 = AD2D_CheckOnBoundary( Current_S_Pos );
%     if(if_on_bdry_0 == 0 && if_on_bdry_1 == 1)
%         Current_S_Mag = 0.5 * Current_S_Mag;
%     end
        
    [~, ~,       FSVec1, FCVec1] = Obtain_System( current_control );
        
    next_sol_S = Current_Sol_S;
    next_sol_C = Current_Sol_C;
    for i=1:n_sub_steps
        FSVec = (n_sub_steps - i) / n_sub_steps * FSVec0 + i / n_sub_steps * FSVec1;

        FCVec0_temp = (1.0 * (n_sub_steps - i) / n_sub_steps + 0.0 )* FCVec0;
        FCVec1_temp = (1.0 * (n_sub_steps - i) / n_sub_steps + 0.0 )* FCVec1;
        FCVec = (n_sub_steps - i) / n_sub_steps * FCVec0_temp + i / n_sub_steps * FCVec1_temp;
        
        KS = M - MatS*dt*theta;
        RhsS = ( MatS*dt*(1-theta) + M )*next_sol_S + FSVec * dt;
        next_sol_S = KS\RhsS;

        KC = M - MatC*dt*theta;
        RhsC = ( MatC*dt*(1-theta) + M )*next_sol_C + FCVec * dt;
        next_sol_C = KC\RhsC;
    end
    Previous_Sol = Current_Sol_S + Current_Sol_C; 
    Current_Sol_S = next_sol_S;
    Current_Sol_C = next_sol_C;
    Current_Sol = Current_Sol_S + Current_Sol_C;
    [reward, scaled_reward] = AD2D_ComputeReward( current_control, Previous_Sol, Current_Sol );

    return;
end


function [MatS, MatC, FSVec, FCVec] = Obtain_System( current_control )
    AD2D_Global_Vars;

    MatS = DS;
    MatC = DC;

    if(Current_S_Pos_Idx <= 0)
        FSVec = zeros(size(FS, 1), 1);
    else
        FSVec = FS(:, Current_S_Pos_Idx);
    end
    %FSVec = ( 0.8 * max(0, 50-Current_Time) / 50 + 0.2 )* FSVec; 
    FSVec = Current_S_Mag * FSVec; 
    
    FCVec = zeros(size(FC, 1), 1);
    if current_control(1) == 1
        MatS = MatS + AL;
        MatC = MatC + AL;
        FCVec = FCVec + FC(:, 1);
    end
    if current_control(2) == 1
        MatS = MatS + AR;
        MatC = MatC + AR;
        FCVec = FCVec + FC(:, 2);
    end
end

