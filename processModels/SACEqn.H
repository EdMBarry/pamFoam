// Scalar transport equation for acetate

fvScalarMatrix SACEqn
(
    fvm::ddt(alphaLiquid, SAC
    + fvm::div(liquidPhase.alphaPhi(), SAC)
    - fvm::laplacian(alphaLiquid*DSND, SAC)
    ==
    - alphaLiquid*rhoAmmonification
    + alphaLiquid*rhoOrganicNitrogenHydrolysis
    - fvm::Sp(dissipationCoeff, SND)
);

SNDEqn.relax();
SNDEqn.solve();
