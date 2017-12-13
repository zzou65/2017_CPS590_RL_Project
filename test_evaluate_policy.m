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

plt_idx = 21:50;
figure(1)
plot(rewards0(plt_idx), '--sk', 'linewidth', 2, 'markersize', 15)
hold on
plot(rewards1(plt_idx), '--ob', 'linewidth', 2, 'markersize', 15)
hold off
set(gca, 'fontsize', 20)
xlabel('Epochs', 'fontsize', 20)
ylabel('Accumulated rewards', 'fontsize', 20)
h=legend('Random', 'Smart');
set(h, 'fontsize', 20)
mean0 = mean(rewards0(plt_idx))
std0 = std(rewards0(plt_idx))
mean1 = mean(rewards1(plt_idx))
std1 = std(rewards1(plt_idx))
   
% Qnet=load('./DQN_FittedQ_Controls/Controler_10000_Samp.mat');
% Qnet=load('./DQN_OnlineQ_Controls/Controler_201_Epoch.mat');
Qnet=load('./DPN_Controls/Controler_1_Epoch.mat');
Qnet=Qnet.net;
Qnet.If_DPN = 1;
Qnet2=load('./DPN_Critic_Controls/Controler_1_Epoch.mat');
Qnet2=Qnet2.net;
Qnet2.If_DPN = 1;

rewards2 = evaluate_policy( Qnet,  r, sample_Path );
rewards3 = evaluate_policy( Qnet2, r, sample_Path );

figure(2)
subplot(1,2,1)
plot(rewards0, '--sk', 'linewidth', 2, 'markersize', 15)
hold on
plot(rewards1, '--ob', 'linewidth', 2, 'markersize', 15)
hold on
plot(rewards2, '--*c', 'linewidth', 2, 'markersize', 15)
hold off
ylim([-5, 1])
subplot(1,2,2)
plot(rewards0, '--sk', 'linewidth', 2, 'markersize', 15)
hold on
plot(rewards1, '--ob', 'linewidth', 2, 'markersize', 15)
hold on
plot(rewards3, '--*c', 'linewidth', 2, 'markersize', 15)
hold off
ylim([-5, 1])



