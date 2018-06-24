fvScalarMatrix SNOEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), SNO)
    + fvm::div(liquidPhase.alphaPhi(), SNO)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, SNO)
    ==
    - alphaLiquid*(1-YH)/2.86/YH*rhoHeterotrophAnoxicGrowth
    + alphaLiquid*rhoAutotrophAerobicGrowth/YA
);
