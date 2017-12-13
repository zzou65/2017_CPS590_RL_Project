
function [prestate, posstate, actionidx, reward, terminal, if_new_play] = DQN_Take_a_Move( Net, epsilon )
    AD2D_Global_Vars;
    if_new_play = 0;
    if(Current_S_Pos_Idx <= 0)
        AD2D_ResetState();
        if_new_play = 1;
        rd = rand(1);
        if(rd > 0.75) % first step is not always taken
            prestate = [];
            posstate = [];
            actionidx = [];
            reward = [];
            terminal = [];
            return;
        end
    end
    if( Net.If_CNN == 0)
        feature = AD2D_Generate_Feature();
        feature = feature(:)';
        feature = normalize(feature, Net.statemu, Net.statesigma);
        prestate = feature;
        Net = nnff(Net, feature, zeros(1, Net.size(end)));
        [~, idx] = max(Net.a{end});
        rd = rand(1);
        if(rd < epsilon)
            actionidx = randi(4);
        else
            actionidx = idx(1);
        end
        current_control = AD2D_Map_ActionIdx_to_Control( actionidx );
        [~, reward]= AD2D_Simulate_one_step(current_control);
        feature = AD2D_Generate_Feature();
        feature = feature(:)';
        feature = normalize(feature, Net.statemu, Net.statesigma);
        posstate = feature;
        if Current_S_Pos_Idx <= 0
            terminal= 1;
        else
            terminal= 0;
        end
    else
        prestate = AD2D_MakeImgFromState( Current_Sol );
        error('Policy not implemented') ;
    end
end

