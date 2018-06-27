// Scalar transport equation for inert solubles

fvScalarMatrix SIEqn
(
    fvm::ddt(alphaLiquid, SI)
    + fvm::div(liquidPhase.alphaPhi(), SI)
    - fvm::laplacian(alphaLiquid*DSI, SI)
    ==
    + alphaLiquid*fSIXS*rhoHydrolysis
    - fvm::Sp(dissipationCoeff, SI)
);

SIEqn.relax();
SIEqn.solve();