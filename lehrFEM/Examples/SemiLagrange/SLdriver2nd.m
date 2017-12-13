cd ../../;
startup;
cd ./Examples/SemiLagrange/;

NREFS=6;
CFL = [0.5];
main_SemiLagrangianSchemes2nd(CFL,NREFS);
NREFS=6;
CFL = [0.7];
main_SemiLagrangianSchemes2nd(CFL,NREFS);
exit
