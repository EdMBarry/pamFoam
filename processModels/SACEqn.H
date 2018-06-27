// Scalar transport equation for acetate

fvScalarMatrix SACEqn
(
    fvm::ddt(alphaLiquid, SAC)
    + fvm::div(liquidPhase.alphaPhi(), SAC)
    - fvm::laplacian(alphaLiquid*DSAC, SAC)
    ==
    - alphaLiquid*rhoAcetateUptake
    + alphaLiquid*fSACXS*rhoHydrolysis
    - fvm::Sp(dissipationCoeff, SAC)
);

SACEqn.relax();
SACEqn.solve();
