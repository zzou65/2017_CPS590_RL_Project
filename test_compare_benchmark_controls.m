clear all
close all
AD2D_Global_Vars;
AD2D_Initialize();
N = 10;
sample_path = AD2D_Simulate_Path( 1, N + 1 );

AD2D_Initialize_PathToFollow( sample_path );
sols1 = zeros(NNode, N);
for i=1:N
    current_control = randi([0, 1], 1, 2);
    [~, ~]= AD2D_Simulate_one_step(current_control);
    sols1(:, i) = Current_Sol;
    if Current_S_Pos_Idx <= 0
        break;
    end
end

AD2D_Initialize_PathToFollow( sample_path );
sols2 = zeros(NNode, N);
for i=1:N
    if Current_S_Pos(1) < 0.5 - 1e-10
        current_control = [1, 0];
    elseif Current_S_Pos(1) > 0.5 + 1e-10
        current_control = [0, 1];
    else
        rd = rand(1);
        if(rd<0.5)
            current_control = [1, 0];
        else
            current_control = [0, 1];
        end
    end
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
    set(gca, 'fontsize', 18)
    subplot(1,2,2)
    plotSurf(MeshObj.Coordinates, sols2(:, k), [-0.05, 0.15], 1)
    set(gca, 'fontsize', 18)
end


