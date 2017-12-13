
function rewards = evaluate_policy( f_Handel, disct_r, samplePath )

    AD2D_Global_Vars;
    [NEpoch, NStep] = size(samplePath);
    rewards = zeros(NEpoch, 1);
   
    for i=1:NEpoch
        fprintf('Epoch %d...\n', i);
        path = samplePath(i, :);
        AD2D_Initialize_PathToFollow( path );
        while( Current_S_Pos_Idx > 0 && Current_Time < NStep - 1 )
            if strcmp( f_Handel, 'random' )
                current_control = randi([0, 1], 1, 2);
            elseif strcmp( f_Handel, 'smart')
                if Current_S_Pos(1) < 0.5 - 1e-10
                    current_control = [1, 0];
                elseif Current_S_Pos(1) > 0.5 + 1e-10
                    current_control = [0, 1];
                else
                    rd = rand(1);
                    if(rd<0.5)
                        current_control = [1, 0];
                    else
                        current_control = [0, 1];
                    end
                    %current_control = [1, 1];
                end
            else
               %error('Policy not implemented') ;
                if( f_Handel.If_CNN == 0)
                    if(f_Handel.If_DPN == 1)
                        idx = DPN_NN_ChooseAction(f_Handel, 4, 1);
                    else
                        feature = AD2D_Generate_Feature();
                        feature = feature(:)';
                        feature = normalize(feature, f_Handel.statemu, f_Handel.statesigma);
                        f_Handel = nnff(f_Handel, feature, zeros(1, f_Handel.size(end)));
                        [~, idx] = max(f_Handel.a{end}(1:4));
                    end
                    current_control = AD2D_Map_ActionIdx_to_Control( idx(1) );
                else
                    prestate = AD2D_MakeImgFromState( Current_Sol );
                    error('Policy not implemented') ;
                end
            end
            % not normalized rewards
            [rewd, ~]= AD2D_Simulate_one_step(current_control);
            %[~, rewd]= AD2D_Simulate_one_step(current_control);
            rewards(i) = rewards(i) + disct_r^Current_Time * rewd;
        end
    end
end
