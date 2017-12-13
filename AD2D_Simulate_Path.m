% sample_Path: NEpoch X MaxNStep

function sample_Path = AD2D_Simulate_Path( NEpoch, MaxNStep )
    AD2D_Global_Vars;
    sample_Path = zeros(NEpoch, MaxNStep);
    for i=1:NEpoch
        path_start_idx = AD2D_AsRdWalk([], SPos, stepSize);
        sample_Path(i, 1) = path_start_idx;
        for j=2:MaxNStep
            new_idx = AD2D_AsRdWalk(sample_Path(i, j-1), SPos, stepSize);
            if(new_idx <= 0)
                break;
            else
                sample_Path(i, j) = new_idx;
            end
        end
    end
end

