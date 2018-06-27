// Scalar transport equation for inorganic phosphorus

fvScalarMatrix SIPEqn
(
    fvm::ddt(alphaLiquid, SIP)
    + fvm::div(liquidPhase.alphaPhi(), SIP)
    - fvm::laplacian(alphaLiquid*SIP, SIP)
    ==
    - fvm::Sp(dissipationCoeff, SIP)
    + alphaLiquid * fSIPXS * rhoHydrolysis
    - alphaLiquid * fPB * YPBPH * rhoAcetateUptake
    - alphaLiquid * fPB * YPBPH * rhoPhotoheterotrophicUptake
    + alphaLiquid * fIPDEC * rhoDecay
);

SIPEqn.relax();
SIPEqn.solve();
