// Scalar transport equation for readily biodegradable substrate

fvScalarMatrix SSEqn
(
    fvm::ddt(alphaLiquid, SS)
    + fvm::div(liquidPhase.alphaPhi(), SS)
    - fvm::laplacian(alphaLiquid*DSS, SS)
    ==
    + alphaLiquid * fSSXS * rhoHydrolysis
    - alphaLiquid * rhoPhotoheterotrophicUptake
    - fvm::Sp(dissipationCoeff, SS)
);

SSEqn.relax();
SSEqn.solve();
