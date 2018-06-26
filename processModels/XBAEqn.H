// Scalar transport equation for active autotrophic biomass

fvScalarMatrix XBAEqn
(
    fvm::ddt(alphaLiquid, XBA)
    + fvm::div(liquidPhase.alphaPhi(), XBA)
    - fvm::laplacian(alphaLiquid*DXBA, XBA)
    ==
      alphaLiquid*rhoAutotrophAerobicGrowth
    - alphaLiquid*rhoAutotrophDecay
    - fvm::Sp(dissipationCoeff, XBA)
);

XBAEqn.relax();
XBAEqn.solve();
