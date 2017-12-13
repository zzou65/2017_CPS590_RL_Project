function Lloc = LOAD_Vol_PWDG_dual(Vertices,ElemData,QuadRule,omega,FHandle,varargin)
%LOAD_VOL_PWDG_DUAL Element load vector for volume data
%
%   LLOC = LOAD_VOL_PWDG_DUAL(VERTICES,ELEMDATA,QUADRULE,OMEGA,FHANDLE)
%   computes the volume contributions of the element load vector for
%   dual (mixed) discontinuous plane waves.
%
%   VERTICES is 3-by-2 or 4-by-2 matrix whose rows contain the vertices of
%   the current element.
%
%   ELEMDATA is a structure that contains at least the fields:
%     NDOFS   The number of degrees of freedom on the current element.
%     DIR     A P-by-2 matrix containing the propagation directions of
%             the plane wave basis functions in its rows.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   The scalar OMEGA is the wave number of the plane waves.
%
%   FHANDLE is a function handle for the load data.  It should take at
%   least the argument x.
%
%   LLOC = LOAD_VOL_PWDG(...,FPARAM1,FPARAM2,...,FPARAMK) also passes the
%   parameters FPARAM1, FPARAM2, ..., FPARAMK to the function handle
%   FHANDLE.
%
%   Example:
%
%     B_Vol = assemLoad_Vol_PDG2_vec(...
%       Mesh,3,@LOAD_Vol_PWDG_dual,qr2,omega,f);
%     B_Imp = assemLoad_Bnd_PDG2_vec(...
%       Mesh,-1,3,@LOAD_Imp_Bnd_PWDG_dual,qr1,omega,gI);
%     B_Dir = assemLoad_Bnd_PDG2_vec(...
%       Mesh,-2,3,@LOAD_Dir_Bnd_PWDG_dual,qr1,omega,gD);
%     B = B_Vol + B_Imp + B_Dir;

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Get load vector from primal method
  Lloc0 = LOAD_Vol_PWDG(Vertices,ElemData,QuadRule,omega,FHandle,varargin);
  
  % Construct load vector for mixed method
  Lloc = [Lloc0/(i*omega);zeros(2*size(Lloc0,1),1)];
  
return