function submemo = AD2D_Subsample_Memory( memo, nsamp, if_add_feature, if_image_format )
%function submemo = AD2D_Subsample_Memory(memo, nsamp, if_add_feature, if_image_format, if_recompute_reward)

    AD2D_Global_Vars;
    if(memo.Image_Format ~= if_image_format)
        error('Image_Format not matching');
    end
    if(memo.Nsamp < nsamp)
        error('Not enough samples in memory');
    end
    if(memo.If_AddFeature == 0 && if_add_feature == 1)
        error('No addfeatures in memory');
    end
    if(memo.If_AddFeature == 0 && if_add_feature == 0)
        st_idx = 1;
    end
    if(memo.If_AddFeature == 1 && if_add_feature == 1)
        st_idx = 1;
    end
    if(memo.If_AddFeature == 1 && if_add_feature == 0)
        st_idx = N_Add_Feature + 1;
    end
    submemo.Nsamp = nsamp;
    submemo.Image_Format = memo.Image_Format;
    submemo.If_AddFeature = if_add_feature;
    submemo.prestate = memo.prestate(1:nsamp, st_idx:end);
    submemo.posstate = memo.posstate(1:nsamp, st_idx:end);
    submemo.rewards = memo.rewards(1:nsamp);
    submemo.actions = memo.actions(1:nsamp);
    submemo.terminal = memo.terminal(1:nsamp);
    submemo.statemu = memo.statemu(1, st_idx:end);
    submemo.statesigma = memo.statesigma(1, st_idx:end);
%{
    if(if_recompute_reward == 1)        
        oldreward = submemo.rewards;
        newreward = oldreward;
        for i=1:nsamp
            prestate = submemo.prestate(i, :) .* submemo.statesigma + submemo.statemu; 
            posstate = submemo.posstate(i, :) .* submemo.statesigma + submemo.statemu;
            control = AD2D_Map_ActionIdx_to_Control( submemo.actions(i) );
            [ ~, newreward(i) ] = AD2D_ComputeReward( control, prestate, posstate);
        end
        submemo.rewards = newreward; 
    end
%}
end

