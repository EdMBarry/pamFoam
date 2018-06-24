fvScalarMatrix XBAEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), XBA)
    + fvm::div(liquidPhase.alphaPhi(), XBA)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, XBA)
    ==
      alphaLiquid*rhoAutotrophAerobicGrowth
    - alphaLiquid*rhoAutotrophDecay
);
