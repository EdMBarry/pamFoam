fvScalarMatrix SOEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), SO)
    + fvm::div(liquidPhase.alphaPhi(), SO)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, SO)
    ==
    - alphaLiquid*(1-YH)/YH*rhoHeterotrophAerobicGrowth
    - alphaLiquid*(4.57-YA)/YA*rhoAutotrophAerobicGrowth
    + alphaLiquid*kLa*(SOsat - SO)
);
