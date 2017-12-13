cd ../../;
startup;
cd ./Examples/SemiLagrange/;

NREFS = 6
CFL = 0.8;
main_SemiLagrangianSchemes(CFL,NREFS)
main_SemiLagrangianCfree(CFL,NREFS)
main_SemiLagrangianPure(CFL,NREFS)
main_SemiLagrangianPureCfree(CFL,NREFS)

CFL = 0.5
main_SemiLagrangianSchemes(CFL,NREFS)
main_SemiLagrangianCfree(CFL,NREFS)
main_SemiLagrangianPure(CFL,NREFS)
main_SemiLagrangianPureCfree(CFL,NREFS)

CFL = 0.25
main_SemiLagrangianSchemes(CFL,NREFS)
main_SemiLagrangianCfree(CFL,NREFS)
main_SemiLagrangianPure(CFL,NREFS)
main_SemiLagrangianPureCfree(CFL,NREFS)
