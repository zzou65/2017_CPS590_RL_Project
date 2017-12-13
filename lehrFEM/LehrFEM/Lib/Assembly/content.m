% ASSEMBLY
%
% Files
%
% Assembly of stiffness matrices
%   assemMat_BFE           - assembles bilinear FE contributions
%
%   assemMat_Bnd_DG        - 
%   assemMat_Bnd_PDG       -
%   assemMat_Bnd_PDG2      -
%   assemMat_Bnd_PDG2_vec  -
%   assemMat_Bnd_Robin     -
%   assemMat_Bnd_hpDG_1D   -
%
%   assemMat_CD_LFV        -
%
%   assemMat_Contr1f       - 
%   assemMat_Contr2f       -
%
%   assemMat_CR            - assembles Crouzeix-Raviart FE contributions
%
%   assemMat_FOSLS_TRP1    -
%
%   assemMat_Inn_DG        -
%   assemMat_Inn_PDG       -
%   assemMat_Inn_PDG2      -
%   assemMat_Inn_PDG2_vec  -
%
%   assemMat_hpDG_1D       -
%
%   assemMat_LFE           - assembles linear FE contributions
%   assemMat_LFE2          - assembles nodal FE contribtutions
%
%   assemMat_LFV           -
%
%   assemMat_Lapl_dual     -
%
%   assemMat_Mass0fD       -
%   assemMat_Mass1fD       -
%   assemMat_Mass2fD       -
%
%   assemMat_P0            - assembles constant FE contributions
%   assemMat_P0_1D         - assembles constant FE contributions (1D)
%
%   assemMat_P1P0          - assembles linear/constant FE contributions
%   assemMat_P1P0_1D       - assembles linear/constant FE contributions (1D)
%
%   assemMat_P1_1D         - assembles linear FE contributions
%
%   assemMat_PBD           -
%
%   assemMat_QFE           - assembles quadratic FE contributions
%                           
%   assemMat_Stokes_CRP0   -
%   assemMat_Stokes_MINIP1 -
%   assemMat_Stokes_P1P0   -
%   assemMat_Stokes_P2P0   -
%   assemMat_Stokes_TH     -
% 
%   assemMat_TopGrad       -
%   assemMat_TopRot        -
%
%   assemMat_UpQFE         -
%
%   assemMat_Vol_DG        -
%   assemMat_Vol_PDG       -
%   assemMat_Vol_PDG2      -
%   assemMat_Vol_PDG2_vec  -
%   assemMat_Vol_hpDG_1D   -
% 
%   assemMat_W1F           - assembles edge elements in 2D
%   assemMat_WRegW1F       - assembles WReg (?) W1F FE contributions
%
%   assemMat_hp            - assembles hp-FEM contributions
%
%   assemMat_minSurf       -
%
%
%
% Assembly of load vectors
%   assemLoad_BFE           - assembles bilinear FE contributions
% 
%   assemLoad_Bnd_DG        -
%   assemLoad_Bnd_PDG       -
%   assemLoad_Bnd_PDG2      -
%   assemLoad_Bnd_PDG2_vec  -
%
%   assemLoad_CR            - assembles Crouzeix-Raviart FE contributions
%
%   assemLoad_FOSLS_TRP1    -
%
%   assemLoad_LFE           - assembles linear FE contributions
%   assemLoad_LFE2          - assembles nodal FE contributions
%   assemLoad_LFE_S1        -
%   assemLoad_LFE_SUPG      -
%
%   assemLoad_LFV           -
%
%   assemLoad_Lapl_dual     -
%
%   assemLoad_P0            - assembles constant FE contributions
%
%   assemLoad_P1_1D         - assembles linear FE contributions (1D)
%
%   assemLoad_PBD           -
%
%   assemLoad_QFE           - assembles quadratic FE contributions
%   assemLoad_QFE_SUPG      -
% 
%   assemLoad_Stokes_CRP0   -
%   assemLoad_Stokes_MINIP1 -
%   assemLoad_Stokes_P1P0   -
%   assemLoad_Stokes_P2P0   -
%   assemLoad_Stokes_TH     -
%
%   assemLoad_Vol_DG        -
%   assemLoad_Vol_PDG       -
%   assemLoad_Vol_PDG2      -
%   assemLoad_Vol_PDG2_vec  -
%
%   assemLoad_W1F           - assembles Whitney 1-forms contributions
%
%   assemLoad_hp            - assembles hp-FEM contributions
%   assemLoad_hpDG_1D       - assembles hpDG-FEM contributions
%
%
%
% Incorporation of Dirichlet boundary conditions
%   assemDir_BFE           - incorporates the Dirichlet boundary conditions into the bilinear FE solution U
%
%   assemDir_CR            - incoporates the Dirichlet boundary conditions into the Crouzeix-Raviart FE solution U
%
%   assemDir_FOSLS_TRP1    -
%
%   assemDir_LFE           - incoporates the Dirichlet boundary conditions into the linear FE solution U
%   assemDir_LFE2          - incoporates the Dirichlet boundary conditions into the (vector-valued) linear FE solution U
%
%   assemDir_LFV           - 
%
%   assemDir_P1_1D         - incoporates the Dirichlet boundary conditions into the linear FE solution U (1D)
%
%   assemDir_PBD           -
%
%   assemDir_QFE           - incoporates the Dirichlet boundary conditions into the quadratic FE solution U
%
%   assemDir_Stokes_CRP0   -
%   assemDir_Stokes_MINIP1 -
%   assemDir_Stokes_P1P0   -
%   assemDir_Stokes_P2P0   -
%   assemDir_Stokes_TH     -
%
%   assemDir_StrRegLFE2    -
%
%   assemDir_W1F           - incoporates the Dirichlet boundary conditions into the edge FE solution U
%
%   assemDir_dual          -
%
%   assemDir_hp            - incoporates the Dirichlet boundary conditions into the hp FE solution U
%   assemDir_hpDG_1D       - incoporates the Dirichlet boundary conditions into the hp FE solution U (1D)
%
%
%
% Incorporation of Neumann boundary conditions
%   assemNeu_BFE   - incoporates the Neumann boundary conditions into the right hand side load vector L for bilinear FE
%   assemNeu_LFE   - incoporates the Neumann boundary conditions into the right hand side load vector L for linear FE
%   assemNeu_LFV   -
%   assemNeu_P1_1D - adds the Neumann boundary conditions at the vertices onto the right-hand side load vector L and the finite element solution U
%   assemNeu_PBD   -
%   assemNeu_QFE   - incoporates the Neumann boundary conditions into the right hand side load vector L for quadratic FE
%   assemNeu_hp    - incoporates the Neumann boundary conditions into the right hand side load vector L for hp FE
%
%
%
% set_
% assemPrec
% assemRKDG
% assemInterp
% assemConv