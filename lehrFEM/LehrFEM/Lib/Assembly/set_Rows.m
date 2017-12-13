function I = set_Rows(ridx,nCols)
% SET_ROWS Row index set.
%
%   I = SET_ROWS(RIDX,NCOLS) row index set for the transformation of an
%   element matrix into an array. RIDX is the set of row indices and NCOLS
%   is the number of columns of the element matrix.
%
%   Example:
%
%   I = set_Rows([1 2 3 4],5);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  I = ridx(ones(nCols,1),:)';
  I = I(:);
  
return