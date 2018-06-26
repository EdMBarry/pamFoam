// Scalar transport equation for inert solubles

fvScalarMatrix SIEqn
(
    fvm::ddt(alphaLiquid, SI)
    + fvm::div(liquidPhase.alphaPhi(), SI)
    - fvm::laplacian(alphaLiquid*SI, SI)
    ==
    + alphaLiquid*fSIXS*rhoHydrolysis
    - fvm::Sp(dissipationCoeff, SNH)
);

SIEqn.relax();
SIEqn.solve();