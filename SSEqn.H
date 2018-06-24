fvScalarMatrix SSEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), SS)
    + fvm::div(liquidPhase.alphaPhi(), SS)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, SS)
    ==
    - alphaLiquid*rhoHeterotrophAerobicGrowth/YH
    - alphaLiquid*rhoHeterotrophAnoxicGrowth/YH
    + alphaLiquid*rhoOrganicHydrolysis
);
