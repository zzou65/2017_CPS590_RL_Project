function x = aux_prec(Fe,Se,P,m,inv_A,GE,T_GE,SN_Handle,varargin)
% AUX_PREC Auxiliary space preconditioner for H(curl)-elliptic problems
%
% Fe -> vector to which preconditioner should be applied
% Se -> edge element Galerkin matrix (for smoother)
% P -> nodal to edge transfer matrix
% GE -> discrete gradient matrix
% T_GE -> transposed of discrete gradient matrix
% SN_Handle -> solver routine in nodal space
%
%   X = AUX_PREC(FE,SE,P,M,SN_HANDLE) computes the value of the auxiliary
%   preconditioner for the matrix Se and the right-hand side Fe.
%
%   Example:
%
%   x = aux_prec(Fe,Se,P,m,SN_Handle);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  
    
    x = zeros(length(Fe),1);
    
    % m steps of symmetric Gauss-Seidel
    
    Le = tril(Se);
    Re = triu(Se);

    for i = 1:m, x = x + Le\(Fe-Se*x); end
    for i = 1:m, x = x + Re\(Fe-Se*x); end
    
    % Auxiliary space correction on Lagrangian finite element space
    x = x + P*SN_Handle(transpose(P)*Fe,varargin{:});
    
    % Auxiliary space correction in discrete potential space
    x = x +  GE*inv_A*T_GE*Fe;
return

