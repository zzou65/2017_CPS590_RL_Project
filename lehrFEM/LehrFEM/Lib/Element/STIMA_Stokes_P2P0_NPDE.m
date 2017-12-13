function Aloc = STIMA_Stokes_P2P0(Vertices,ElemInfo,nu,QuadRule,varargin)
% Computation of element matrix for P2-P0 finite element discretization of 2D Stokes problem
% \texttt{Vertices} passes the location of the vertices of the triangle
% \texttt{nu} is the viscosity parameter
% \texttt{QuadRule} specifies local quadrature rule, see Rem.~\ref{rem:lfemquad}
% The function returns a \Blue{$14\times14$} dense matrix
Aloc = zeros(14,14); % Preallocate memory
% Compute element mapping
bK = Vertices(1,:); BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
inv_BK_t = transpose(inv(BK)); det_BK = abs(det(BK));
% Compute gradients element shape functions and their values at quadrature points
grad_N = grad_shap_QFE(QuadRule.x);
grad_N(:,1:2) = grad_N(:,1:2)*inv_BK_t;
grad_N(:,3:4) = grad_N(:,3:4)*inv_BK_t;
grad_N(:,5:6) = grad_N(:,5:6)*inv_BK_t;
grad_N(:,7:8) = grad_N(:,7:8)*inv_BK_t;
grad_N(:,9:10) = grad_N(:,9:10)*inv_BK_t;
grad_N(:,11:12) = grad_N(:,11:12)*inv_BK_t;
% The first 6 rows/columns of the element matrix correspond to the \Blue{$x_1$}-component of the
% velocity. the corresponding block of the element matrix agrees with that for \Blue{$-\Delta$}
% discretized by means of quadratic Lagrangian finite elements. The local shape functions are
% described in Ex.~\ref{ex:quadLFE}.
Aloc(1,1) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,1:2),2))*det_BK;
Aloc(1,2) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,3:4),2))*det_BK;
Aloc(1,3) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,5:6),2))*det_BK;
Aloc(1,4) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,7:8),2))*det_BK;
Aloc(1,5) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,9:10),2))*det_BK;
Aloc(1,6) = nu*sum(QuadRule.w.*sum(grad_N(:,1:2).*grad_N(:,11:12),2))*det_BK;
Aloc(2,2) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,3:4),2))*det_BK;
Aloc(2,3) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,5:6),2))*det_BK;
Aloc(2,4) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,7:8),2))*det_BK;
Aloc(2,5) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,9:10),2))*det_BK;
Aloc(2,6) = nu*sum(QuadRule.w.*sum(grad_N(:,3:4).*grad_N(:,11:12),2))*det_BK;
Aloc(3,3) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,5:6),2))*det_BK;
Aloc(3,4) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,7:8),2))*det_BK;
Aloc(3,5) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,9:10),2))*det_BK;
Aloc(3,6) = nu*sum(QuadRule.w.*sum(grad_N(:,5:6).*grad_N(:,11:12),2))*det_BK;
Aloc(4,4) = nu*sum(QuadRule.w.*sum(grad_N(:,7:8).*grad_N(:,7:8),2))*det_BK;
Aloc(4,5) = nu*sum(QuadRule.w.*sum(grad_N(:,7:8).*grad_N(:,9:10),2))*det_BK;
Aloc(4,6) = nu*sum(QuadRule.w.*sum(grad_N(:,7:8).*grad_N(:,11:12),2))*det_BK;
Aloc(5,5) = nu*sum(QuadRule.w.*sum(grad_N(:,9:10).*grad_N(:,9:10),2))*det_BK;
Aloc(5,6) = nu*sum(QuadRule.w.*sum(grad_N(:,9:10).*grad_N(:,11:12),2))*det_BK;
Aloc(6,6) = nu*sum(QuadRule.w.*sum(grad_N(:,11:12).*grad_N(:,11:12),2))*det_BK;
% the same for the \Blue{$x_2$}-component of the velocity
Aloc(7,7) = Aloc(1,1);  Aloc(7,8) = Aloc(1,2);  Aloc(7,9) = Aloc(1,3);
Aloc(7,10) = Aloc(1,4);  Aloc(7,11) = Aloc(1,5);  Aloc(7,12) = Aloc(1,6);
Aloc(8,8) = Aloc(2,2);  Aloc(8,9) = Aloc(2,3);  Aloc(8,10) = Aloc(2,4);
Aloc(8,11) = Aloc(2,5);  Aloc(8,12) = Aloc(2,6);  Aloc(9,9) = Aloc(3,3);
Aloc(9,10) = Aloc(3,4);  Aloc(9,11) = Aloc(3,5);  Aloc(9,12) = Aloc(3,6);
Aloc(10,10) = Aloc(4,4);  Aloc(10,11) = Aloc(4,5);  Aloc(10,12) = Aloc(4,6);
Aloc(11,11) = Aloc(5,5);  Aloc(11,12) = Aloc(5,6);  Aloc(12,12) = Aloc(6,6);
% Interaction of pressure shape functon (constant \Blue{$\equiv 1$}) with velocity: evaluation of
% local bilinear form \Blue{$\blf{b}_K$}. 
% First for \Blue{$x_1$}-components, then for 
Aloc(1,13) = sum(QuadRule.w.*grad_N(:,1))*det_BK;
Aloc(2,13) = sum(QuadRule.w.*grad_N(:,3))*det_BK;
Aloc(3,13) = sum(QuadRule.w.*grad_N(:,5))*det_BK;
Aloc(4,13) = sum(QuadRule.w.*grad_N(:,7))*det_BK;
Aloc(5,13) = sum(QuadRule.w.*grad_N(:,9))*det_BK;
Aloc(6,13) = sum(QuadRule.w.*grad_N(:,11))*det_BK;
% Next for \Blue{$x_2$}-components of velocity
Aloc(7,13) = sum(QuadRule.w.*grad_N(:,2))*det_BK;
Aloc(8,13) = sum(QuadRule.w.*grad_N(:,4))*det_BK;
Aloc(9,13) = sum(QuadRule.w.*grad_N(:,6))*det_BK;
Aloc(10,13) = sum(QuadRule.w.*grad_N(:,8))*det_BK;
Aloc(11,13) = sum(QuadRule.w.*grad_N(:,10))*det_BK;
Aloc(12,13) = sum(QuadRule.w.*grad_N(:,12))*det_BK;
% Entry corresponding to zero mean multiplier
Aloc(13,14) = det_BK;
% Fill in lower triangular part
tri = triu(Aloc); Aloc = tri+tril(tri',-1);
return