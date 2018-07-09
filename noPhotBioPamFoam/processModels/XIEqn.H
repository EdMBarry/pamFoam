// Scalar transport equation for particulate products arising from biomass decay

fvScalarMatrix XIEqn
(
    fvm::ddt(alphaLiquid, XI)
    + fvm::div(liquidPhase.alphaPhi(), XI)
    - fvm::laplacian(alphaLiquid*DXI, XI)
    ==
    - fvm::Sp(dissipationCoeff, XI)
    + alphaLiquid*fXIXS*rhoHydrolysis
);

XIEqn.relax();
XIEqn.solve();