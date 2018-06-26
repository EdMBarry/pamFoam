// Define concentrations to be used in limiting rates
dimensionedScalar zeroRate = dimensionedScalar("zeroRate", dimensionSet(1,-3,-1,0,0,0,0), 0);
dimensionedScalar smallConc = dimensionedScalar("smallConc", dimensionSet(1,-3,0,0,0,0,0), SMALL);

// Define rate limiting, competition, and inhibiting expressions
volScalarField inhibFA = max(kFA / (SIN * kFA), smallConc); 
volScalarField limitingSIN = max(SIN / (kSIN * SIN), smallConc);
volScalarField limitingSIP = max(SIP / (kSIP * SIP), smallConc);
volScalarField compSAC = max(SAC / (SS + SAC), smallConc);
volScalarField compSS = max(SS / (SAC + SS), smallConc);
volScalarField limitingG = max(GLambda_0 / (kG + GLambda_0), smallConc);

// Define multiplier to zero reaction rates in pure gas regions
volScalarField mult = 1 - pos(alphaGas - 0.5);

// Define process rates
rhoHeterotrophAerobicGrowth = mult * max(muH * SS/(KS + SS) * SO/(KOH + SO) * XBH, zeroRate);
rhoHeterotrophAnoxicGrowth = mult * max(muH * SS/(KS + SS) * KOH/(KOH + SO) * SNO/(KNO + SNO) * etag * XBH, zeroRate);
rhoPhotoheterotrophicGrowth = mult * max(muPH * XPB * IFA * IIN * IIP * iG * IAC * ICA, zeroRate);
rhoDecay = mult * max(kd * XPB, zeroRate);
rhoHydrolysis = mult * max(kh * XS, zeroRate);