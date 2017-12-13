% Runs experiment 2 for a given set of parameters, then
% saves all numerical data to exp2_lfv.mat.


    exp2_epsrange = [1e-0 5e-1 1e-1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4];
    exp2_stdmesh = [];
    exp2_solutions_1 = [];
    refs = 4;
    
    refrange = 0:6;
    exp2_meshes = [];
    exp2_solutions_2 = [];
    exp2_errors_L1 = [];
    exp2_errors_L2 = [];
    exp2_errors_Linf = [];
    exp2_hrange = [];

    for e=1:size(exp2_epsrange,2)
        disp(['Running test ' num2str(e) ' of ' num2str(size(exp2_epsrange,2)+size(refrange,2))]);
        [errs, mw, msh, u] = exp_2(exp2_epsrange(e), refs, 0, 0);
        exp2_stdmesh = msh;
        exp2_solutions_1(e).sol = u;
    end



    for r=1:size(refrange,2)
        disp(['Running test ' num2str(size(exp2_epsrange,2)+r) ' of ' num2str(size(exp2_epsrange,2)+size(refrange,2))]);
        [errs, mw, msh, u] = exp_2(0, refrange(r), 0, 0);
        exp2_errors_L1(r) = errs(1);
        exp2_errors_L2(r) = errs(2);
        exp2_errors_Linf(r) = errs(3);
        exp2_solutions_2(r).sol = u;

        exp2_meshes(r).Coordinates = msh.Coordinates;
        exp2_meshes(r).Elements = msh.Elements;
        exp2_meshes(r).ElemFlag = msh.ElemFlag;
        exp2_meshes(r).Edges = msh.Edges;
        exp2_meshes(r).Vert2Edge = msh.Vert2Edge;
        exp2_meshes(r).Max_Nodes = msh.Max_Nodes;
        exp2_meshes(r).BdFlags = msh.BdFlags;
        exp2_meshes(r).CenterPoints = msh.CenterPoints;
        exp2_meshes(r).Type = msh.Type;
        exp2_meshes(r).MidPoints = msh.MidPoints;

        exp2_hrange(r) = mw;
    end

    

    save 'exp2_lfv.mat' exp2_*
