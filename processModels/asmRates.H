// Define concentrations to be used in limiting rates
dimensionedScalar zeroRate = dimensionedScalar("zeroRate", dimensionSet(1,-3,-1,0,0,0,0), 0);
dimensionedScalar smallConc = dimensionedScalar("smallConc", dimensionSet(1,-3,0,0,0,0,0), SMALL);

// Define multiplier to zero reaction rates in pure gas regions
volScalarField mult = 1 - pos(alphaGas - 0.5);

// Define process rates
rhoHeterotrophAerobicGrowth = mult * max(muH * SS/(KS + SS) * SO/(KOH + SO) * XBH, zeroRate);
rhoHeterotrophAnoxicGrowth = mult * max(muH * SS/(KS + SS) * KOH/(KOH + SO) * SNO/(KNO + SNO) * etag * XBH, zeroRate);
rhoAutotrophAerobicGrowth = mult * max(muA * SNH/(KNH + SNH) * SO/(KOA + SO) * XBA, zeroRate);
rhoHeterotrophDecay = mult * max(bH * XBH, zeroRate);
rhoAutotrophDecay = mult * max(bA * XBA, zeroRate);
rhoAmmonification = mult * max(ka * SND * XBH, zeroRate);
rhoOrganicHydrolysis = mult * max(kh * XS/(KX*XBH + max(XS, smallConc)) * (SO/(KOH + SO) + etah * KOH/(KOH + SO) * SNO/(KNO + SNO)) * XBH, zeroRate);
rhoOrganicNitrogenHydrolysis = mult * max(rhoOrganicHydrolysis * XND / max(XS, smallConc), zeroRate);
