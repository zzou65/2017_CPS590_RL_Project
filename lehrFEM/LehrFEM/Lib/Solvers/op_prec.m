function x = op_prec(b,it,cyc,m,MG_Data)
% OP_PREC Operator preconditioner.
%
%   X = OP_PREC(B,IT,CYC,M,MG_DATA) computes the value of the operator 
%   preconditioner given by the matrices MG_DATA.A and the right hand
%   side B using a multigrid solver for fast evaluation. 
%
%   IT specifies the number of multigrid iterations to be performed for the
%   preconditioner.
%
%   CYC denotes the type of cycle to be used at each multigrid iteration:
%    1 V-cycles with M1 pre- and M2 post-smoothing steps.
%    2 W-cycles with M1 pre- and M2 post-smoothing steps.
%
%   The struct MG_DATA should at least contain the following fields:
%    A is a 1-by-LVL cell array containing the preconditioner matrix on
%      each level.
%    P 1-by-(LVL-1) cell array containing the prolongation matrix on each
%      level.
%
%   Example:
%
%   x = op_prec(b,1,1,MG_Data);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Run multigrid solver
 
  x = b;
  lvl = size(MG_Data.A,2);
  for i = 1:it
    x = mg_step(x,b,lvl,cyc,m,MG_Data);
  end
    
return
      
%%% Multigrid steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

function x = mg_step(x,b,lvl,cyc,m,MG_Data)
% MG Multigrid solver.
%
%   X = MG_STEP(X0,B,LVL,CYC,M,MG_DATA) performs one multigrid step
%   from the initial guess X0, the load vector B and the stiffness
%   matrices given by cell array MG_DATA.A using the multigrid method.
%
%   LVL denotes the mesh level for which the multigrid step is supposed to
%   be run.
%
%   CYC denotes the type of cycle to be used for multigrid method:
%    1 V-cycles with M1 pre- and M2 post-smoothing steps.
%    2 W-cycles with M1 pre- and M2 post-smoothing steps.
%
%   The struct MG_DATA should at least contain the following fields:
%    A is a 1-by-LVL cell array containing the stiffness matrix on each
%      level.
%    P 1-by-(LVL-1) cell array containing the prolongation matrix on each
%      level.
%
%   Example:
%
%   x = mg(x,b,lvl,cyc,m,MG_Data);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  if(lvl == 1)
    
    % Use direct solver on the coarse mesh  
    
    if(~isempty(MG_Data.A{lvl}))
      x = MG_Data.A{lvl}\b;
    end
    
  else
           
    % Extract matrices
      
    A = MG_Data.A{lvl};
    P = MG_Data.P{lvl-1};
    L = tril(A);
    R = triu(A);   
    
    % Relax M times (Gauss-Seidel smoother)  
      
    for i = 1:m
      x = x+L\(b-A*x);  
    end
    
    % Correct fine grid solution (V or W cycles)
    
    res = b-A*x;    
    res = transpose(P)*res;
    cor = zeros(size(res));
    for i = 1:cyc
      cor = mg_step(cor,res,lvl-1,cyc,m,MG_Data);
    end    
    x = x+P*cor;
    
    % Relax M times (Gauss-Seidel smoother)
    
    for i = 1:m
      x = x+R\(b-A*x);  
    end
    
  end

return