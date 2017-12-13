% memo is initial memory, ONLY used for normalization
% output is n action + 1 value
function net = DPN_Train_NN ( archi, opts, memo, if_use_baseline, n_epoch, n_trajectory, disct_r, max_traj_len)
    AD2D_Global_Vars;

    if(memo.Image_Format == 1)
        error('NN needs nonimage input')
    end
    if(if_use_baseline == 1 && archi(end) ~= 5)
        error('Wrong output archi')
    end
    if(if_use_baseline == 0 && archi(end) ~= 4)
        error('Wrong output archi')
    end
        
    % setup network
    net = nnsetup( archi );
    net.output = 'linear';
    net.learningRate = opts.learningRate;
    net.momentum = opts.momentum;
    net.dropoutFraction = opts.dropoutFraction;
    net.weightPenaltyL2 = opts.weightPenaltyL2;
    net.If_CNN = 0;
    net.statemu = memo.statemu;
    net.statesigma = memo.statesigma;
    net.If_DPN = 1;
    net.If_BaseLine = if_use_baseline;
    
    % setup exploration strategy
    init_epsilon = opts.init_epsilon;
    deca_epsilon = opts.deca_epsilon;

    for i=1:n_epoch
        fprintf('Training epoch %d ...\n', i);
        saveloc = sprintf('./DPN_Controls/Controler_%d_Epoch.mat', i);
        current_epsilon = max(0.05, init_epsilon - deca_epsilon * (i-1));
        trajbatch = cell(1, n_trajectory);
        for j=1:n_trajectory
            fprintf('traj: %d ...\n', j);
            trajmemo = Simulate_A_Trajectory(net, disct_r, max_traj_len, current_epsilon);
            trajbatch{j} = trajmemo;
        end
       [batchsize, trainX, trainY] = process_memo_into_trdata(trajbatch, disct_r, if_use_baseline);
        tr_opts.numepochs =  1;
        tr_opts.batchsize = batchsize;
        net = nntrain(net, trainX, trainY, tr_opts);
        if(mod(i, 50) == 1)
            save(saveloc, 'net');
        end
    end
end

function [bsize, trainX, trainY] = process_memo_into_trdata(trajbatch, disct_r, if_use_b)
    n_traj = length(trajbatch);
    bsize = 0; 
    trainX = [];
    trainY = [];
    n_out = size(trajbatch{1}.nnout, 2);
    if(if_use_b == 1)
        n_act = n_out - 1;
    else
        n_act = n_out;
        %all_reward = [];
        %for i=1:n_traj
        %    all_reward = [all_reward; trajbatch{i}.disct_rewards(:)];
        %end
        %all_reward_mean = mean(all_reward);
        %all_reward_std  = std (all_reward);
    end
    for i=1:n_traj
        memo = trajbatch{i};
        memoX = memo.prestate;
        memoY = zeros(memo.nstep, n_out);
        for j=1:memo.nstep
            actidx = memo.actions(j);
            disct_rewd = memo.disct_rewards(j);
            if(if_use_b == 1)
                advantage = disct_rewd - memo.nnout(j, n_out);
            else
                %advantage = ( disct_rewd - all_reward_mean) / all_reward_std;
                advantage = disct_rewd;
            end
            disct_fac = disct_r^(j-1);
            for k=1:n_act
                if(k==actidx)
                    delta = disct_fac * advantage * (1.0 - memo.actionprobs(j, k));
                else
                    delta = disct_fac * advantage * (0.0 - memo.actionprobs(j, k));
                end
                memoY(j, k) = memo.nnout(j, k) + delta;
            end
            if(if_use_b == 1)
                memoY(j, n_out) = disct_rewd;
            end
        end
        bsize = bsize + memo.nstep;
        trainX = [trainX; memoX];
        trainY = [trainY; memoY];
    end
end

function trajmemo = Simulate_A_Trajectory(nn, disct_r, maxN, cur_ep)
    AD2D_Global_Vars;
    AD2D_ResetState();
    trajmemo.prestate=[];
    trajmemo.posstate=[];
    trajmemo.nnout=[];
    trajmemo.actionprobs=[];
    trajmemo.actions=[];
    trajmemo.rewards=[];
    trajmemo.disct_rewards=[];
    trajmemo.terminal=[];
    trajmemo.nstep = 0;
    while(Current_S_Pos_Idx > 0 && trajmemo.nstep < maxN)
        trajmemo.nstep = trajmemo.nstep + 1;
        feature = AD2D_Generate_Feature();
        feature = feature(:)';
        feature = normalize(feature, nn.statemu, nn.statesigma);
        prestate = feature;
        trajmemo.prestate = [trajmemo.prestate; prestate];
        nn = nnff(nn, feature, zeros(1, nn.size(end)));
        trajmemo.nnout=[trajmemo.nnout; nn.a{end}(:)'];
        [actionprobs, actionidx] = chooseaction(nn.a{end}, 4, cur_ep);
        trajmemo.actionprobs=[trajmemo.actionprobs; actionprobs(:)'];
        trajmemo.actions=[trajmemo.actions; actionidx];
        current_control = AD2D_Map_ActionIdx_to_Control( actionidx );
        [~, rewd]= AD2D_Simulate_one_step(current_control);
        trajmemo.rewards=[trajmemo.rewards; rewd];
        feature = AD2D_Generate_Feature();
        feature = feature(:)';
        feature = normalize(feature, nn.statemu, nn.statesigma);
        posstate = feature;
        trajmemo.posstate = [trajmemo.posstate; posstate];
        if Current_S_Pos_Idx <= 0
            terminal= 1;
        else
            terminal= 0;
        end
        trajmemo.terminal = [trajmemo.terminal; terminal];
    end
    trajmemo.disct_rewards = compute_disct_rewards(trajmemo.rewards, disct_r);
end
        
function disct_rewards = compute_disct_rewards( rewards, disct_r )
    G = rewards;
    for i = length(rewards)-1:-1:1
        G(i) = G(i) + disct_r * G(i+1); 
    end
    %%% should I normalize the G here?
    %G = (G - mean(G))./std(G);
    disct_rewards = G;
end
 
function [actionprobs, actionidx] = chooseaction(nnoutvec, naction, cur_ep)
        expout = exp(nnoutvec(1:naction));
        actionprobs = expout./sum(expout);
        probsinterval = cumsum(actionprobs);
        probsinterval = [0; probsinterval(:)];
        if(abs(probsinterval(end)-1) > 1e-10)
            error('Check me: total action prob not 1');
        end
        rd = rand(1);
        if(rd < cur_ep)
            actionidx = randi(4);
        else
            rd = rand(1);
            actionidx = 0;
            for i=1:naction
                if(rd > probsinterval(i) && rd < probsinterval(i+1))
                    actionidx = i;
                    break;
                end
            end
            if(actionidx == 0)
                [~, actionidx] = max(actionprobs);
            end
        end
end
    

