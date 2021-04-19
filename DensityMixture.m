function [ density ] = DensityMixture(mfAS,mfOrg,mfWater,rhoAS,rhoOrg,kk)

    rhoSolInt = 1.2; % initial guess value of solution density

    for i=1:10

        mfBYrhoSum = 1./rhoSolInt - mfWater./0.998; % calculated from initialized density
        mfBYrhoSumActual = mfAS./rhoAS + mfOrg./rhoOrg;

        rhoSol = 1./(mfAS./rhoAS + mfOrg./rhoOrg + mfWater./0.998); % assuming ideal behavior in dilute solutions

        error(i,kk) = (mfBYrhoSum-mfBYrhoSumActual).^2./mfBYrhoSum.^2;

        rhoSolInt = rhoSol;

    end

    assignin('base', 'ErrorDensity', error);
  
    density = rhoSol;
    
end

