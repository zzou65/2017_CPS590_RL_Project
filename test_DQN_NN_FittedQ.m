clear all
close all
AD2D_Global_Vars;
AD2D_Initialize();
addpath(genpath('./DeepLearnToolbox'));

Ntr = 100000;
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
nsamp = 100000;
memo = AD2D_Subsample_Memory(memo, nsamp, 1, 0);

saveloc = sprintf('./DQN_FittedQ_Controls/Controler_%d_Samp.mat', nsamp);
archi = [size(memo.prestate, 2), 50, 4];
opts.learningRate = 0.02;
opts.momentum = 0.2;
opts.dropoutFraction = 0.0;
opts.weightPenaltyL2 = 1e-4;
n_epoch = 2000;
batchsize = 2500;
disct_r = 0.95;
Qnet = DQN_Train_NN_FittedQ ( archi, opts, memo, n_epoch, batchsize, disct_r);
save(saveloc, 'Qnet')

