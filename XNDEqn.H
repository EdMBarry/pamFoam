fvScalarMatrix XNDEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), XND)
    + fvm::div(liquidPhase.alphaPhi(), XND)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, XND)
    ==
      alphaLiquid*(iXB-fP*iXP)*rhoHeterotrophDecay
    + alphaLiquid*(iXB-fP*iXP)*rhoAutotrophDecay
    - alphaLiquid*rhoOrganicNitrogenHydrolysis
);
