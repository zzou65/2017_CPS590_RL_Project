function J = set_Cols(cidx,nRows)
% SET_COLS Column index set.
%
%   J = SET_COLS(CIDX,NROWS) column index set for the transformation of an
%   element matrix into an array. CIDX is the set of column indices and
%   NROWS the number of rows of the element matrix.
%
%   Example:
%
%   I = set_Cols([1 2 3 4],5);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  J = cidx(ones(1,nRows),:);
  J = J(:);
  
return