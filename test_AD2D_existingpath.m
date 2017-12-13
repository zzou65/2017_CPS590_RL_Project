close all

N = 20;
sample_Path = AD2D_Simulate_Path( 1, N + 1 );
AD2D_Initialize_PathToFollow( sample_path );

sols = zeros(NNode, N);
reward_hist = [];
FS_mag_hist = [];
pos_hist = Current_S_Pos;
control_hist = [];
for i=1:N
    current_control = randi([0, 1], 1, 2);
    [rewd]= AD2D_Simulate_one_step(current_control);
    sols(:, i) = Current_Sol;
    reward_hist = [reward_hist, rewd];
    FS_mag_hist = [FS_mag_hist, Current_S_Mag];
    control_hist = [control_hist; current_control];
    figure(i)
    plotSurf(MeshObj.Coordinates, Current_Sol, [], 1)
    if Current_S_Pos_Idx > 0
        pos_hist = [pos_hist; SPos(Current_S_Pos_Idx, :)];
    else
        break;
    end
end

figure(100)
plot(pos_hist(:, 1), pos_hist(:, 2), '--o')
xlim([0, 1])
ylim([0, 1])

figure(101)
subplot(1,3,1)
plot(control_hist(:, 1), 'ok', 'markersize', 15)
hold on
plot(control_hist(:, 2), 'xr', 'markersize', 15)
hold off
subplot(1,3,2)
plot(reward_hist, '--or')
subplot(1,3,3)
plot(FS_mag_hist, '--or')

