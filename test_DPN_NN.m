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

if_use_baseline = 0;
archi = [size(memo.prestate, 2), 50, 4];
opts.learningRate = 0.05;
opts.momentum = 0.1;
opts.dropoutFraction = 0.0;
opts.weightPenaltyL2 = 1e-4;
opts.init_epsilon = 1.0;
opts.deca_epsilon = 0.03;
n_epoch = 1001;
n_trajectory = 10;
disct_r = 0.95;
max_traj_len = 30;

net = DPN_Train_NN ( archi, opts, memo, if_use_baseline, n_epoch, n_trajectory, disct_r, max_traj_len);
save('DPN_NNControl.mat', 'net');

