// Define concentrations to be used in limiting rates
dimensionedScalar zeroRate = dimensionedScalar
                            (
                                "zeroRate",
                                dimensionSet(1,-3,-1,0,0,0,0),
                                VSMALL
                            );

dimensionedScalar smallConc = dimensionedScalar
                            (
                                "smallConc",
                                dimensionSet(1,-3,0,0,0,0,0),
                                SMALL
                            );

dimensionedScalar zeroRateAU = dimensionedScalar
                            (
                                "zeroRateAU",
                                dimensionSet(0,-3,-1,0,1,0,0),
                                VSMALL
                            );

// Define rate limiting, competition, and inhibiting expressions
Info << "reading inhibFA" << nl << endl;
inhibFA = max((KFA / (SIN + KFA)), SMALL);

Info << "reading limitingSAC" << nl << endl;
limitingSAC = max(SAC / (KSAC + SAC), SMALL);

Info << "reading limitingSS" << nl << endl;
limitingSS = max(SS / (KSS + SS), SMALL);

limitingSIC = max(SIC / (KSIC + SIC), SMALL);

limitingSH2 = max(SH2 / (KSH2 + SH2), SMALL);

Info << "reading limitingSIN" << nl << endl;
limitingSIN = max(SIN / (KSIN + SIN), SMALL);

Info << "reading limitingSIP" << nl << endl;
limitingSIP = max(SIP / (KSIP + SIP), SMALL);

Info << "reading limitingG" << nl << endl;
limitingG = max(G850/ (KG + G850), SMALL);

Info << "reading compSAC" << nl << endl;
volScalarField compSAC = max((SS + smallConc)/ (SS + SAC + smallConc), SMALL);

Info << "reading compSS" << nl << endl;
volScalarField compSS = max((SAC + smallConc)/ (SAC + SS + smallConc), SMALL);

// Define multiplier to zero reaction rates in pure gas regions
Info << "reading mult" << nl << endl;
volScalarField mult = 1 - pos(alphaGas - 0.9);

// Define process rates

Info << "reading rhoAcetateUptake" << nl << endl;
rhoAcetateUptake = mult
                * max(muAC*XPB*inhibFA*limitingSIN*limitingSIP
                * limitingG*compSS*limitingSAC,
                zeroRate
                );


Info << "reading rhoHetUptake" << nl << endl;
rhoPhotoheterotrophicUptake = mult * max(muPH * XPB * inhibFA * limitingSIN
                                * limitingSIP * limitingG * limitingSS
                                * compSAC, zeroRate);

Info << "Reading rhoChemoheterotrophicUptake" << nl << endl;
rhoChemoheterotrophicUptake = mult * max(muCH * XPB * inhibFA * limitingSIN
                                * limitingSIP * limitingSS, zeroRate);


Info << "Reading rhoAutotrophicUptake" << nl << endl;
rhoAutotrophicUptake = mult * max(muAU * XPB * inhibFA * limitingSIN
                            * limitingSIP * limitingG * limitingSIC
                            * limitingSH2, zeroRateAU);

Info << "reading rhoDecay" << nl << endl;
rhoDecay = mult * max(kd * XPB, zeroRate);

Info << "reading rhoHyd" << nl << endl;
rhoHydrolysis = mult * max(kh * XS, zeroRate);
