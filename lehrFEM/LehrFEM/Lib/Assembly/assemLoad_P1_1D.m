function L = assemLoad_P1_1D(Coordinates,QuadRule,FHandle,varargin)
% ASSEMLOAD_P1_1D Assemble linear FE contributions.
%
%   L = ASSEMLOAD_P1_1D(COORDINATES,QUADRULE,FHANDLE) assembles the global
%   load vector for the load data given by the function handle FHANDLE.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_P1_1D(COORDINATES,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function
%   handle FHANDLE.
%
%   Example:
%
%   Coordinates = transpose(0:.01:1);
%   FHandle = @(x,varargin)8*ones(size(x));
%   L = assemLoad_P1_1D(Coordinates,gauleg(0,1,5),FHandle);
%
%   See also shap_1D.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,1);
  
  % Preallocate memory
  
  Lloc = zeros(2,1);
  L = zeros(nCoordinates,1);
  
  % Precompute shape function values
  
  N = shap_P1_1D(QuadRule.x);
  
  for i = 1:(nCoordinates-1)

    % Compute element mapping
      
    h = abs(Coordinates(i+1)-Coordinates(i));
    x = Coordinates(i) + h*QuadRule.x;
    
    % Evaluate load data at quadrature points
    
    FVal = FHandle(x,varargin{:});
    
    % Add contributions to global load vector
    
    L(i) = L(i) + sum(QuadRule.w.*FVal.*N(:,1))*h;
    L(i+1) = L(i+1) + sum(QuadRule.w.*FVal.*N(:,2))*h;
    
  end

return