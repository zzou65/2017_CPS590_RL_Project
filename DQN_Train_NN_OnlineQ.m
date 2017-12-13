% memo is initial memory
function net = DQN_Train_NN_OnlineQ ( archi, opts, memo, n_epoch, n_steps_each_epoch, batchsize, disct_r)
    if(memo.Image_Format == 1)
        error('NN needs nonimage input')
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

    % setup exploration strategy
    init_epsilon = opts.init_epsilon;
    deca_epsilon = opts.deca_epsilon;

    oldnn = net;
    tr_opts.numepochs =  1;
    tr_opts.batchsize = batchsize;

    memosize = memo.Nsamp;
    n_total = 0;
    n_play = 1;
    for i=1:n_epoch
        fprintf('Training epoch %d ...\n', i);
        saveloc = sprintf('./DQN_OnlineQ_Controls/Controler_%d_Epoch.mat', i);
        current_epsilon = max(0.05, init_epsilon - deca_epsilon * (i-1));
        for j=1:n_steps_each_epoch
            [prestate, posstate, actionidx, reward, terminal, if_new_play] = DQN_Take_a_Move( oldnn, current_epsilon );
            if(isempty(prestate))
                continue;
            end
            if(if_new_play > 1e-10)
                n_play = n_play + 1;
            end
            n_total = n_total + 1;
            pos = mod(n_total, memosize);
            if(pos == 0)
                pos = memosize;
            end
            % unpdate the memory 
            memo.prestate(pos, :) = prestate;
            memo.posstate(pos, :) = posstate;
            memo.actions(pos) = actionidx;
            memo.rewards(pos) = reward;
            memo.terminal(pos) = terminal;
            
            batch_idx = randperm(memosize);
            batch_idx = batch_idx(1:batchsize);
            submemo = sample_from_memo(memo, batch_idx);
            %target = NN_Generate_Y( oldnn, submemo, disct_r );
            target = NN_Generate_Y( oldnn, net, submemo, disct_r );
            net = nntrain(net, submemo.prestate, target, tr_opts);
        end
        oldnn = net;
        net.n_play = n_play;
        if(mod(i, 100) == 1)
            save(saveloc, 'net');
        end
    end
end

function submemo = sample_from_memo(memo, idx)
    submemo.Nsamp = length(idx);
    submemo.prestate = memo.prestate(idx, :);
    submemo.posstate = memo.posstate(idx, :);
    submemo.rewards = memo.rewards(idx);
    submemo.actions = memo.actions(idx);
    submemo.terminal = memo.terminal(idx);
end

function targ = NN_Generate_Y( oldnet, net, memo, r )
       net = nnff(   net, memo.prestate, zeros(memo.Nsamp,    net.size(end)));
    targ =  net.a{end};
    oldnet = nnff(oldnet, memo.posstate, zeros(memo.Nsamp, oldnet.size(end)));
    yy = oldnet.a{end};
    yymax = max(yy, [], 2);
    for i = 1:memo.Nsamp
        reward = memo.rewards(i);
        % only change one action, other Q-values stay the same
        if memo.terminal(i) == 1
            targ(i, memo.actions(i)) = reward;
        else
            targ(i, memo.actions(i)) = reward + r * yymax(i);
        end
    end
end


