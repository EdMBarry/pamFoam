fvScalarMatrix XSEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), XS)
    + fvm::div(liquidPhase.alphaPhi(), XS)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, XS)
    ==
      alphaLiquid*(1-fP)*rhoHeterotrophDecay
    + alphaLiquid*(1-fP)*rhoAutotrophDecay
    - alphaLiquid*rhoOrganicHydrolysis
);
