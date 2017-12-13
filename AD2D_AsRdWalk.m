% starting from [1/4, 1/2]
% idx = -1 if out of bounds
% grid always [0, 1] X [0, 1]
% u/d/l/r + 0.1 with p 0.35, 0.15, 0.15, 0.35

function next_idx = AD2D_AsRdWalk(current_idx, grid, sSize)
    if(isempty(current_idx))
        next_pos = [1/4, 1/2];
        next_idx = check_idx(next_pos(1), next_pos(2), grid);
    elseif(current_idx <= 0)
        next_idx = -1;
    else
        next_pos = grid(current_idx, :);
        r = rand(1);
        if r < 0.35
            next_pos(2) = next_pos(2) + sSize;
        elseif r < 0.70
            next_pos(1) = next_pos(1) + sSize;
        elseif r < 0.85
            next_pos(2) = next_pos(2) - sSize;
        else
            next_pos(1) = next_pos(1) - sSize;
        end
        next_idx = check_idx(next_pos(1), next_pos(2), grid);
    end 
end

function idx = check_idx(x, y, grid)
    if_out = check_boundary(x, y, [0, 1], [0, 1]);
    if(if_out)
        idx = -1;
        return;
    end
    
    [grid_type, dx, dy] = check_grid_type(grid);
    nx = round(1/dx) + 1;
    ny = round(1/dy) + 1;
    m = round(x/dx);
    n = round(y/dy);
    if(grid_type == 1)
        idx = m * ny + n + 1;
    else
        idx = n * nx + m + 1;
    end
end

function if_out = check_boundary(x, y, xrange, yrange)
    if_out = 0;
    if ( x < xrange(1) - 1e-10 || x > xrange(2) + 1e-10 )
        if_out = 1;
    end
    if ( y < yrange(1) - 1e-10 || y > yrange(2) + 1e-10 )
        if_out = 1;
    end
    return;
end

function [grid_type, dx, dy] = check_grid_type(grid)
    if abs(grid(2,1)-grid(1,1)) < 1e-10
        grid_type = 1;
        dy = abs(grid(2,2)-grid(1,2));
        x0 = grid(1, 1);
        x1 = grid(2, 1);
        i = 3;
        while abs(x1 - x0) < 1e-10
            x0 = x1;
            x1 = grid(i, 1);
            i = i + 1;
        end
        dx = abs(x1 - x0);
    else
        grid_type = 2;
        dx = abs(grid(2,1)-grid(1,1));
        y0 = grid(1, 2);
        y1 = grid(2, 2);
        i = 3;
        while abs(y1 - y0) < 1e-10
            y0 = y1;
            y1 = grid(i, 2);
            i = i + 1;
        end
        dy = abs(y1 - y0);
    end
end
