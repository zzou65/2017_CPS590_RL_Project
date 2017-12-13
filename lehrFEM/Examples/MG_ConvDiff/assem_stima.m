function A_fd = assem_stima(mesh,c,v)
%ASSEM_STIMA assemble stiffness matrix for convection-diffusion
%   
%   A = ASSEM_STIMA(MESH,C,V) assembles a stiffness matrix for the
%   operator
%     -C*\Delta + V*\nabla
%   on the mesh MESH using an upwind quadrature scheme.  C can be either a
%   scalar or a function handle; V must be a function handle.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  if(isa(c,'numeric'))
    A = c*assemMat_LFE(mesh,@STIMA_Lapl_LFE);
  else
    A = assemMat_LFE(mesh,@STIMA_Heat_LFE,P1O2(),c);
  end
  M = assemMat_MassZeroD(mesh);
  D_fd = assemMat_LFE(mesh,@LIEUP_FD,v);

  A_fd = A+M*D_fd;
  
return