% Runs experiment 1 for a given set of parameters, then
% saves all numerical data to exp1_lfv.mat.


    exp1_epsrange = [1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1];
    refrange = 1:6;
    exp1_meshes = [];
    exp1_solutions = [];
    exp1_errors_L1 = [];
    exp1_errors_L2 = [];
    exp1_errors_Linf = [];
    exp1_hrange = [];

    for r=1:size(refrange,2)
        for e=1:size(exp1_epsrange,2)
            disp(['Running test ' num2str((e-1)+size(exp1_epsrange,2)*(r-1)+1) ' of ' num2str(size(exp1_epsrange,2)*size(refrange,2))]);
            [errs, mw, msh, u] = exp_1(exp1_epsrange(e), refrange(r), 0, 0);
            exp1_errors_L1(e,r) = errs(1);
            exp1_errors_L2(e,r) = errs(2);
            exp1_errors_Linf(e,r) = errs(3);
            exp1_solutions(e,r).sol = u;
        end

        exp1_meshes(r).Coordinates = msh.Coordinates;
        exp1_meshes(r).Elements = msh.Elements;
        exp1_meshes(r).ElemFlag = msh.ElemFlag;
        exp1_meshes(r).Edges = msh.Edges;
        exp1_meshes(r).Vert2Edge = msh.Vert2Edge;
        exp1_meshes(r).Max_Nodes = msh.Max_Nodes;
        exp1_meshes(r).BdFlags = msh.BdFlags;
        exp1_meshes(r).CenterPoints = msh.CenterPoints;
        exp1_meshes(r).Type = msh.Type;
        exp1_meshes(r).MidPoints = msh.MidPoints;
        exp1_hrange(r) = mw;
    end

    

    save 'exp1_lfv.mat' exp1_*
