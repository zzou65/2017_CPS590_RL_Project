
function image_proj_mat = AD2D_MakeImgProjMat(n_pixels_per_d)
    AD2D_Global_Vars;
    npd = n_pixels_per_d;
    dx = 1/npd;
    dy = 1/npd;
    image_proj_mat = zeros(npd^2, NNode);

    NodesDivFloors = zeros(NNode, 2);
    for i=1:NNode
        NodesDivFloors(i, 1) = floor( MeshObj.Coordinates(i, 1) / dx );
        NodesDivFloors(i, 2) = floor( MeshObj.Coordinates(i, 2) / dy );
        NodesDivFloors(i, 1) = round( max( min(NodesDivFloors(i, 1), npd - 1), 0 ) );
        NodesDivFloors(i, 2) = round( max( min(NodesDivFloors(i, 2), npd - 1), 0 ) );
        pixel_pos = NodesDivFloors(i, 1) * npd + NodesDivFloors(i, 2) + 1;
        image_proj_mat(pixel_pos, i) = 1.0;
    end
    for i=1:size(image_proj_mat, 1)
        image_proj_mat(i, :) = image_proj_mat(i, :)/sum(image_proj_mat(i,:));
    end

