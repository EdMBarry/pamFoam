// Scalar transport equation for slowly biodegradable substrate

fvScalarMatrix XSEqn
(
    fvm::ddt(alphaLiquid, XS)
    + fvm::div(liquidPhase.alphaPhi(), XS)
    - fvm::laplacian(alphaLiquid*DXS, XS)
    ==
    - fvm::Sp(dissipationCoeff, XS)
    + alphaLiquid*rhoDecay
    - alphaLiquid*rhoHydrolysis
);

XSEqn.relax();
XSEqn.solve();
