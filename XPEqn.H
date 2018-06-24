fvScalarMatrix XPEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), XP)
    + fvm::div(liquidPhase.alphaPhi(), XP)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, XP)
    ==
      alphaLiquid*fP*rhoHeterotrophDecay
    + alphaLiquid*fP*rhoAutotrophDecay
);
