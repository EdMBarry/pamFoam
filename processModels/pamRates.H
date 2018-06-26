// Define concentrations to be used in limiting rates
dimensionedScalar zeroRate = dimensionedScalar
                            (
                                "zeroRate",
                                dimensionSet(1,-3,-1,0,0,0,0),
                                0
                            );

dimensionedScalar smallConc = dimensionedScalar
                            (
                                "smallConc",
                                dimensionSet(1,-3,0,0,0,0,0),
                                SMALL
                            );

// Define rate limiting, competition, and inhibiting expressions
volScalarField inhibFA = max(kFA / (SIN * kFA), smallConc); 

volScalarField limitingSAC = max(SAC / (KSAC + SAC), smallConc);
volScalarField limitingSS = max(SS / (KSS + SS), smallConc);
volScalarField limitingSIN = max(SIN / (kSIN + SIN), smallConc);
volScalarField limitingSIP = max(SIP / (kSIP + SIP), smallConc);
volScalarField limitingG = max(GLambda_0 / (kG + GLambda_0), smallConc);

volScalarField compSAC = max(SS / (SS + SAC), smallConc);
volScalarField compSS = max(SAC / (SAC + SS), smallConc);

// Define multiplier to zero reaction rates in pure gas regions
volScalarField mult = 1 - pos(alphaGas - 0.5);

// Define process rates
rhoAcetateUptake = mult
                * max
                (
                    muAC * XPB * inhibFA * limitingSIN * limitingSIP
                    * limitingG * compSS * limitingSAC,
                    zeroRate
                );


rhoPhotoheterotrophicUptake = mult * max(muPH * XPB * inhibFA * limitingSIN
                                * limitingSIP * limitingG * limitingSS 
                                * compSAC, zeroRate);
rhoDecay = mult * max(kd * XPB, zeroRate);
rhoHydrolysis = mult * max(kh * XS, zeroRate);