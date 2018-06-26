// Scalar transport equation for slowly biodegradable substrate

fvScalarMatrix XSEqn
(
    fvm::ddt(alphaLiquid, XS)
    + fvm::div(liquidPhase.alphaPhi(), XS)
    - fvm::laplacian(alphaLiquid*DXS, XS)
    ==
      alphaLiquid*(1-fP)*rhoHeterotrophDecay
    + alphaLiquid*(1-fP)*rhoAutotrophDecay
    - alphaLiquid*rhoOrganicHydrolysis
    - fvm::Sp(dissipationCoeff, XS)
);

XSEqn.relax();
XSEqn.solve();
