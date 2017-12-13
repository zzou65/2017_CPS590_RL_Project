function x = permute_smoother(x,A,b,smoother,varargin)
%PERMUTE_SMOOTHER reorder indices before applying smoother
%
%   X = PERMUTE_SMOOTHER(X,A,B,SMOOTHER,PER1,PER2,...,PERN)
%   applies the smoother SMOOTHER to iteratively solve A*X=B with each of
%   the permutations PER1,PER2,...,PERN.  A permutation PERK must be a
%   2-by-M matrix if M is the dimension of A.  The first row, PERK(1,:), is
%   the actual permutation of the indices applied before the smoother and
%   PERK(2,:) is the inverse permutation applied afterwards.  The inverse
%   permutation INVPER can be calculated from the forward permutation PER
%   by [dummy,INVPER] = sort(PER).

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  n = length(varargin);

  for k=1:n

    % define permutation
    
    per = varargin{k}(1,:);
    iper = varargin{k}(2,:);
    % [dummy,iper] = sort(per);
    
    % apply smoother
    
    xp = smoother(x(per),A(per,per),b(per));
    
    % apply inverse permutation
    
    x = xp(iper);

  end
  
return