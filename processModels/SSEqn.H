// Scalar transport equation for readily biodegradable substrate

fvScalarMatrix SSEqn
(
    fvm::ddt(alphaLiquid, SS)
    + fvm::div(liquidPhase.alphaPhi(), SS)
    - fvm::laplacian(alphaLiquid*DSS, SS)
    ==
    - fvm::Sp(dissipationCoeff, SS)
    + alphaLiquid * fSSXS * rhoHydrolysis
    - alphaLiquid * rhoPhotoheterotrophicUptake
    - alphaLiquid * rhoChemoheterotrophicUptake
);

SSEqn.relax();
SSEqn.solve();
