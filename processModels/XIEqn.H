// Scalar transport equation for particulate products arising from biomass decay

fvScalarMatrix XIEqn
(
    fvm::ddt(alphaLiquid, XI)
    + fvm::div(liquidPhase.alphaPhi(), XI)
    - fvm::laplacian(alphaLiquid*XI, XI)
    ==
    + alphaLiquid*fXIXS*rhoHydrolysis
    - fvm::Sp(dissipationCoeff, XI)
);

XIEqn.relax();
XIEqn.solve();