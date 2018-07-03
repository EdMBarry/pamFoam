// Scalar transport equation for photoheterotrophic biomass

fvScalarMatrix XPBEqn
(
    fvm::ddt(alphaLiquid, XPB)
    + fvm::div(liquidPhase.alphaPhi(), XPB)
    - fvm::laplacian(alphaLiquid*DXPB, XPB)
    ==
    - fvm::Sp(dissipationCoeff, XPB)
    - alphaLiquid*rhoDecay
    + alphaLiquid*YPBPH*rhoAcetateUptake
    + alphaLiquid*YPBPH*rhoPhotoheterotrophicUptake
    + alphaLiquid*YPBCH*rhoChemoheterotrophicUptake
    + alphaLiquid*YPBAU*rhoAutotrophicUptake
);

XPBEqn.relax();
XPBEqn.solve();
