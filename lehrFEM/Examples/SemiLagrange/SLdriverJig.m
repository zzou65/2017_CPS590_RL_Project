cd ../../;
startup;
cd ./Examples/SemiLagrange/;

NREFS=6;
CFL = [0.25];
main_SemiLagrangianSchemesJig(CFL,NREFS);
NREFS=6;
CFL = [0.5];
main_SemiLagrangianSchemesJig(CFL,NREFS);
NREFS=6;
CFL = [0.8];
main_SemiLagrangianSchemesJig(CFL,NREFS);
exit
