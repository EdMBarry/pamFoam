// Scalar transport equation for inorganic phosphorus

fvScalarMatrix SIPEqn
(
    fvm::ddt(alphaLiquid, SIP)
    + fvm::div(liquidPhase.alphaPhi(), SIP)
    - fvm::laplacian(alphaLiquid*SIP, SIP)
    ==
    alphaLiquid * fIPXS * rhoHydrolysis
    - alphaLiquid * fPB * YPBPH * rhoAcetateUptake
    - alphaLiquid * fPB * YPBPH * rhoPhotoheterotrophicUptake
    + alphaLiquid * fIPDEC * rhoDecay
    - fvm::Sp(dissipationCoeff, SIP)
);

SIPEqn.relax();
SIPEqn.solve();
