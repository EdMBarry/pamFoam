// Scalar transport equation for slowly biodegradable organic nitrogen

fvScalarMatrix SNDEqn
(
    fvm::ddt(alphaLiquid, SND)
    + fvm::div(liquidPhase.alphaPhi(), SND)
    - fvm::laplacian(alphaLiquid*DSND, SND)
    ==
    - alphaLiquid*rhoAmmonification
    + alphaLiquid*rhoOrganicNitrogenHydrolysis
    - fvm::Sp(dissipationCoeff, SND)
);

SNDEqn.relax();
SNDEqn.solve();
