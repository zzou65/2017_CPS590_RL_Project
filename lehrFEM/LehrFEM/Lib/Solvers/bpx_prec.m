function x = bpx_prec(b,ML_Data,varargin)
% BPX_PREC Bramble-Pasciak-Xu preconditioner.
%
%   X = BPX_PREC(B,ML_DATA) computes the value of the Bramble-Pasciak-Xu
%   preconditioner for the matrix specified by its diagonal ML_DATA.D and
%   the right-hand side B.
%
%   The struct ML_DATA should at least contain the following fields:
%    D is a 1-by-LVL cell array containing the diagonal of the stiffness
%      matrix on each level.
%    P 1-by-(LVL-1) cell array containing the prolongation matrix on each
%      level.
%
%   Example:
%
%   x = bpx_prec(b,ML_Data);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  
 
  x = bpx_step(b,size(ML_Data.D,2),ML_Data);

return

%%% BPX step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = bpx_step(x,lvl,ML_Data)
% BPX_STEP Recursive step of Bramble-Pasciak-Xu preconditioner.
%
%   X = BPX_STEP(X,LVL,ML_DATA) recursivly computes the value of the
%   Bramble-Pasciak-Xu preconditioner for the matrix specified by its
%   diagonal L_DATA.D and the right-hand side X on the level LVL.   

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  D = ML_Data.D{lvl};
  if(lvl > 1)
    P = ML_Data.P{lvl-1};   
    x = x./D + P*bpx_step(transpose(P)*x,lvl-1,ML_Data);
  else
    if(~isempty(D))
      x = x./D;   
    end
  end

return