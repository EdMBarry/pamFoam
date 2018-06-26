// Scalar transport equation for photoheterotrophic biomass

fvScalarMatrix XPBEqn
(
    fvm::ddt(alphaLiquid, XPB)
    + fvm::div(liquidPhase.alphaPhi(), XPB)
    - fvm::laplacian(alphaLiquid*DXPB, XPB)
    ==
     alphaLiquid*YPBPH*rhoAcetateUptake
     alphaLiquid*YPBPH*rhoPhotoheterotrophicUptake
    - alphaLiquid*rhoDecay
    - fvm::Sp(dissipationCoeff, XPB)
);

XPBEqn.relax();
XPBEqn.solve();
