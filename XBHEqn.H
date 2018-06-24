fvScalarMatrix XBHEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), XBH)
    + fvm::div(liquidPhase.alphaPhi(), XBH)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, XBH)
    ==
      alphaLiquid*rhoHeterotrophAerobicGrowth
    + alphaLiquid*rhoHeterotrophAnoxicGrowth
    - alphaLiquid*rhoHeterotrophDecay
);
