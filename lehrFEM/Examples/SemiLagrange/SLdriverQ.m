cd ../../;
startup;
cd ./Examples/SemiLagrange/;

NREFS=6;
CFL = [0.8,0.7,0.6,0.5,0.4,0.3,0.3];
main_SemiLagrangianSchemesTime(CFL,NREFS);
exit
