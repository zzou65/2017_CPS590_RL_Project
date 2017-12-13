
function net = DQN_Train_NN_FittedQ ( archi, opts, memo, n_epoch, batchsize, disct_r)
    if(memo.Image_Format == 1)
        error('NN needs nonimage input')
    end
    net = nnsetup( archi );
    net.output = 'linear';
    net.learningRate = opts.learningRate;
    net.momentum = opts.momentum;
    net.dropoutFraction = opts.dropoutFraction;
    net.weightPenaltyL2 = opts.weightPenaltyL2;
    net.If_CNN = 0;
    net.statemu = memo.statemu;
    net.statesigma = memo.statesigma;

    oldnn = net;
    tr_opts.numepochs =  1;
    tr_opts.batchsize = batchsize;
    for i=1:n_epoch
        fprintf('Training epoch %d ...\n', i);
        target = NN_Generate_Y( oldnn, memo, disct_r );
        net = nntrain(oldnn, memo.prestate, target, tr_opts);
        oldnn = net;
    end
end

function targ = NN_Generate_Y( oldnet, memo, r )
    oldnet = nnff(oldnet, memo.prestate, zeros(memo.Nsamp, oldnet.size(end)));
    targ = oldnet.a{end};
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

