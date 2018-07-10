// Scalar transport equation for inert solubles

fvScalarMatrix SIEqn
(
    fvm::ddt(alphaLiquid, SI)
    + fvm::div(liquidPhase.alphaPhi(), SI)
    - fvm::laplacian(alphaLiquid*DSI, SI)
    ==
    - fvm::Sp(dissipationCoeff, SI)
    + alphaLiquid*fSIXS*rhoHydrolysis
);

SIEqn.relax();
SIEqn.solve();