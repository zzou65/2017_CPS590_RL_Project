% control [0, 0], [1, 0], [0, 1], [1, 1]

function [ reward, scaled_reward ] = AD2D_ComputeReward( current_control, current_sol, next_sol )
    AD2D_Global_Vars;
    potential_0  = potential_func( current_sol );
    potential_1  = potential_func( next_sol );

%{    
    reward = 2 * potential_1;
    if( potential_1 >= -1e-10 )
        reward = reward + 0.10;
    elseif( potential_1 > 0.5 * potential_0 )
        reward = reward + 0.03;
    elseif( potential_1 < 1.5 * potential_0 )
        reward = reward  -0.03;
    else
        reward = reward + 0.00;
    end
%}

    reward = potential_1;
    reward = 2 * reward + 1 * (0.95 * potential_1 - potential_0);
    if( potential_1 >= -1e-10 )
        reward = reward + 0.10;
    end
    
    
    if current_control(1)*current_control(2) < 1-1e-10
        reward = reward - current_control(1) * costL;  
        reward = reward - current_control(2) * costR; 
    else
        reward = reward - current_control(1) * costL * 1.5;  
        reward = reward - current_control(2) * costR * 1.5;
    end

   scaled_reward = min( max(reward, reward_low), reward_hig ); 
   scaled_reward = ( scaled_reward - reward_low ) / (reward_hig - reward_low);
end


function potential = potential_func( sol )
    AD2D_Global_Vars;
    err_vec = zeros(NNode, 1);
    for i = 1:NNode
        val = sol(i);
        if(val >= cost_thre_low && val <= cost_thre_hig)
            err_vec(i) = 0.0;
        elseif(val > cost_thre_hig)
            err_vec(i) = val - cost_thre_hig;
        else
            err_vec(i) = cost_thre_low - val;
        end
    end
    potential = - sqrt(err_vec' * M * err_vec);
end

