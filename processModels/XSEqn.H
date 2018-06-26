// Scalar transport equation for slowly biodegradable substrate

fvScalarMatrix XSEqn
(
    fvm::ddt(alphaLiquid, XS)
    + fvm::div(liquidPhase.alphaPhi(), XS)
    - fvm::laplacian(alphaLiquid*DXS, XS)
    ==
      alphaLiquid*rhoDecay
    - alphaLiquid*rhoHydrolysis
    - fvm::Sp(dissipationCoeff, XS)
);

XSEqn.relax();
XSEqn.solve();
