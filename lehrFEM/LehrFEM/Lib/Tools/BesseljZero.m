function root = BesseljZero(N,X0)
% BESSELJZERO the zero point of the besselj function
% 
% ROOT = BESSELJZERO(NUM,X0) returns the zero point of the zero point of the
% Nth order of besselj function around X0.
%
% Example
%
% root = BesseljZero(0,1);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Initialize constant
    
    FHandle = @(x)besselj(N,x);
    nPoints = length(X0);
    
    % Preallocate memory
    
    root = zeros(nPoints,1);
    
    % Search the zero point
    
    for i = 1:nPoints
        root(i) = fzero(FHandle,X0(i));
    end

return
