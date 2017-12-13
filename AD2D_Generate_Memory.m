% image_format nXnXnsamp, [0, 1], nsamp last dim
% non_image_format nsampXn^2, normalized, nsamp first dim

function Memo = AD2D_Generate_Memory( Ntr, maxNepoch, maxNstep, if_image_format )
    AD2D_Global_Vars;
    n_samp = 0;
    nepoch = 0;
    if(if_image_format == 1)
        Memo.Image_Format = 1;
        Memo.If_AddFeature = 0;
        Memo.Nsamp = Ntr;
        Memo.prestate = zeros(N_Pixels_Per_Dim, N_Pixels_Per_Dim, Ntr);
        Memo.posstate = zeros(N_Pixels_Per_Dim, N_Pixels_Per_Dim, Ntr);
        Memo.actions = zeros(1, Ntr);
        Memo.rewards = zeros(1, Ntr);
        Memo.terminal = zeros(1, Ntr);
    else
        Memo.Image_Format = 0;
        Memo.If_AddFeature = If_Use_AddFeature;
        Memo.Nsamp = Ntr;
        Memo.prestate = zeros(Ntr, N_Feature);
        Memo.posstate = zeros(Ntr, N_Feature);
        Memo.actions = zeros(Ntr, 1);
        Memo.rewards = zeros(Ntr, 1);
        Memo.terminal = zeros(Ntr, 1);
    end
        
    while n_samp < Ntr && nepoch < maxNepoch
        AD2D_ResetState();
        nepoch = nepoch + 1;
        for i=1:maxNstep
            prestate = AD2D_MakeImgFromState( Current_Sol );
            prefeatu = AD2D_Generate_Feature(); 
            action = randi([0, 1], 1, 2);
            [~, s_rewd]= AD2D_Simulate_one_step( action);
            posstate = AD2D_MakeImgFromState( Current_Sol );
            posfeatu = AD2D_Generate_Feature();
            action = AD2D_Map_Control_to_ActionIdx( action );
            rd = rand(1);
            if( i > 1 || (i == 1 && rd > 0.75) ) % first step is not always stored
                n_samp = n_samp + 1;
                if(if_image_format == 1)
                    Memo.prestate(:, :, n_samp) = prestate;
                    Memo.posstate(:, :, n_samp) = posstate;
                else
                    Memo.prestate(n_samp, :) = prefeatu(:)';
                    Memo.posstate(n_samp, :) = posfeatu(:)';
                end
                Memo.actions(n_samp) = action;
                Memo.rewards(n_samp) = s_rewd;
            end
            
            if i == maxNstep
                fprintf('Reached simu step bdd ... breaking ...n_samp: %d ...\n', n_samp)
            end
            if Current_S_Pos_Idx <= 0
                fprintf('Reached terminal state ... breaking ...n_samp: %d ...\n', n_samp)
                Memo.terminal(n_samp) = 1;
                break;
            end
            if n_samp >= Ntr
                fprintf('Memory full ... breaking ...\n')
                break;
            end
        end
    end
    if(if_image_format ~= 1)
        % normalize if not image format
        [Memo.prestate, Memo.statemu, Memo.statesigma] = zscore(Memo.prestate);
        Memo.posstate = normalize(Memo.posstate, Memo.statemu, Memo.statesigma);
    end
end

