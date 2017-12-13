clear all
close all

AD2D_Global_Vars;
AD2D_Initialize();
addpath(genpath('./DeepLearnToolbox'));

r = 0.95;
NEpoch = 200;
MaxNStep = 30;
str = sprintf('test_path_rd_smart_res_%d_runs.mat', NEpoch);
if(exist(str, 'file'))
    load(str);
else
    sample_Path = AD2D_Simulate_Path( NEpoch, MaxNStep + 1);
    save('test_sample_path.mat', 'sample_Path');
    rewards0 = evaluate_policy( 'random', r, sample_Path );
    rewards1 = evaluate_policy( 'smart', r, sample_Path );
    save(str, 'sample_Path', 'rewards0', 'rewards1');
end

mean0 = mean(rewards0);
stdd0 = std(rewards0);
mean1 = mean(rewards1);
stdd1 = std(rewards1);

%QnetNRange = [5000, 10000, 20000, 50000, 100000];
QnetNRange = [1, 201, 401, 601, 801, 1001];
mean_plot2 = [];
stdd_plot2 = [];
for i=1:length(QnetNRange)
    nsamp = QnetNRange(i);
    %saveloc = sprintf('./DQN_FittedQ_Controls/Controler_%d_Samp.mat', nsamp);
    %saveloc = sprintf('./DQN_OnlineQ_Controls/Controler_%d_Epoch.mat', nsamp);
    saveloc = sprintf('./DPN_Controls/Controler_%d_Epoch.mat', nsamp);
    Qnet=load(saveloc);
    %Qnet=Qnet.Qnet;
    Qnet=Qnet.net;
    %Qnet.If_DPN = 0;
    Qnet.If_DPN = 1;
    rewards2 = evaluate_policy( Qnet, r, sample_Path );
    mean_plot2 =[mean_plot2; mean(rewards2)];
    stdd_plot2 =[stdd_plot2;  std(rewards2)];
end

%save('./DQN_FittedQ_Controls/Controler_Rewards.mat', 'QnetNRange', 'mean_plot2', 'stdd_plot2');
%save('./DQN_OnlineQ_Controls/Controler_Rewards.mat', 'QnetNRange', 'mean_plot2', 'stdd_plot2');
save('./DPN_Controls/Controler_Rewards.mat', 'QnetNRange', 'mean_plot2', 'stdd_plot2');

%%%%%%
mean_plot0 = mean0 * ones(length(QnetNRange), 1);
stdd_plot0 = stdd0 * ones(length(QnetNRange), 1);
mean_plot1 = mean1 * ones(length(QnetNRange), 1);
stdd_plot1 = stdd1 * ones(length(QnetNRange), 1);

%mean_plot2(end) = mean_plot2(end-1) * 1.05;
figure(1)
errorbar(QnetNRange, mean_plot0, stdd_plot0, '-k', 'linewidth', 2)
hold on
errorbar(QnetNRange, mean_plot1, stdd_plot1, '-b', 'linewidth', 2)
hold on
errorbar(QnetNRange, mean_plot2, stdd_plot2, '-c', 'linewidth', 2)
hold on
errorbar(QnetNRange, dpn_ac.mean_plot2, dpn_ac.stdd_plot2, '-r', 'linewidth', 2)
hold off
xlim([-10, 1010])
ylim([-1.5, 1])
set(gca, 'fontsize', 20)
xlabel('N_{epoch}', 'fontsize', 20)
title('Accumulated rewards: mean and std', 'fontsize', 20, 'fontweight', 'normal')
h=legend('Random', 'Smart', 'DPN', 'DPN-AC');
set(h, 'fontsize', 18);


