clear all
close all
AD2D_Global_Vars;
AD2D_Initialize();
addpath(genpath('./DeepLearnToolbox'));

Ntr = 5000;
if_image_format = 0;
str = sprintf('memo_Ntr_%d_IfImage_%d.mat', Ntr, if_image_format);
maxNepoch = 100000;
maxNstep = 30;
if(exist(str, 'file'))
    load(str)
else
    memo = AD2D_Generate_Memory( Ntr, maxNepoch, maxNstep, if_image_format );
    save(str, 'memo');
end

archi = [size(memo.prestate, 2), 50, 4];
opts.learningRate = 0.01;
opts.momentum = 0.2;
opts.dropoutFraction = 0.0;
opts.weightPenaltyL2 = 1e-4;
opts.init_epsilon = 1.0;
opts.deca_epsilon = 0.01;
n_epoch = 1001;
n_steps_each_epoch = 50;
batchsize = 200;
disct_r = 0.95;
Qnet = DQN_Train_NN_OnlineQ ( archi, opts, memo, n_epoch, n_steps_each_epoch, batchsize, disct_r);
save('NNControl_Online.mat', 'Qnet')
