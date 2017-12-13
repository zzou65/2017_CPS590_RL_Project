clear all
close all
AD2D_Global_Vars;
AD2D_Initialize();
N = 20;
sample_path = AD2D_Simulate_Path( 1, N + 1 );
AD2D_Initialize_PathToFollow( sample_path );
addpath(genpath('./DeepLearnToolbox'));
Qnet=load('./DQN_NN_Controlers_OnlineQ/NNControler_OnlineQ_351_Epoch.mat');
Qnet=Qnet.net;

sols1 = zeros(NNode, N);
for i=1:N
    [~, ~, ~, ~, ~, ~] = DQN_Take_a_Move( Qnet, 0.0 );
    sols1(:, i) = Current_Sol;
    if Current_S_Pos_Idx <= 0
        break;
    end
end

AD2D_Initialize_PathToFollow( sample_path );
sols2 = zeros(NNode, N);
for i=1:N
    current_control = randi([0, 1], 1, 2);
    [~, ~]= AD2D_Simulate_one_step(current_control);
    sols2(:, i) = Current_Sol;
    if Current_S_Pos_Idx <= 0
        break;
    end
end

for k=1:i
    figure(k)
    subplot(1,2,1)
    plotSurf(MeshObj.Coordinates, sols1(:, k), [-0.05, 0.15], 1)
    subplot(1,2,2)
    plotSurf(MeshObj.Coordinates, sols2(:, k), [-0.05, 0.15], 1)
end


