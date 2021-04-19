clc
clear all

starttimer = tic;

now = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm');

global STw Vw R T
STw = 0.073; % N/m
Vw = 18e-6; % m^3/mol
R = 8.314; % J/mol-K
T = 298.15; % K

prompt = 'ENTER details for size distribution --> [1,#particles,geomean of dia (um),std dev] OR [2,#particles,Min dia (um),Max dia (um)] : ';
[InputSizeDistribution] = input(prompt);
k=InputSizeDistribution(1);n=InputSizeDistribution(2); mu=InputSizeDistribution(3); sigma=InputSizeDistribution(4);

prompt = 'Salt? [ENTER 1 for Ammonium sulfate OR 2 for Sodium chloride]: ';
salt = input(prompt);

prompt = 'Organic? [ENTER 1 for Glutaric acid OR 2 for Succinic acid OR 3 for sucrose]: ';
org = input(prompt);

prompt = 'Black carbon present? [ENTER y or n]: ';
core = input(prompt,'s');

prompt = 'If y, ENTER 1/2/3 for coated sphere, else ENTER 0 for homogeneous sphere: ';
opt = input(prompt);

prompt = '[ENTER 1 for constant BC core dia OR 2 for constant ratio Dshell/Dcore OR 3 for constant shell thickness OR 0 for no BC core]: ';
BC_core = input(prompt);

prompt = 'Enter number of angles for S1 and S2 phase amplitudes in range from 0 to pi/2 (anywhere between 10 and 1000; e.g. 91): ';

angleNo = input(prompt);
interval = 90/(angleNo-1);
theta = 0:interval:180; % size should be (2*angleNo-1); the additive (interval) is calculated by (90/(angleNo-1))

% lambdaRange = 350:50:750; % RI available for these wavelengths only in this visible spectrum
lambdaRange = 550;
NN = size(lambdaRange');

for kkkk = 1:NN
    
    lambda(1,1,1,kkkk) = lambdaRange(kkkk);
    
    %     OIRrange = [1e-100,0.01,0.1,0.2,0.5,1,2,5,10,100];
%         OIRrange = [1e-100,0.15,0.35,0.55,0.75,0.95];
%         OIRrange = 1e-100;
    
    fr_org = 0.05:0.1:0.95;
%     fr_org = [0.25,0.55,0.95];
%     fr_org = 1e-100;
    
    OIRrange = fr_org./(1-fr_org);
    
    OIRrange = [1e-100,OIRrange,1e+100];
    
    NNN = size(OIRrange');
    
    for kkk = 1:NNN
        
        OIR(1,1,kkk,kkkk) = OIRrange(kkk);
                        RHrange = 85:1:99;
        %                         RHrange = 99:-1:85;
%                                 RHrange = [85,95,99];
%         RHrange = 99;
        N = size(RHrange');
        
%         rhoAS = 1.8; % Density of ammonium sulfate in ug/nl
%         rhoAS = 2.2; % Density of sodium chloride in ug/nl
        rhoSalt = [1.8;2.2];
        
        rhoOrg = [1.424;1.572;1.6];
        
%             MolWtAS = 132.14; % Molecular weight of ammonium sulfate (g/mol)
%             MolWtAS = 58.44; % Molecular weight of sodium chloride (g/mol)
        MolWtSalt = [132.14;58.44];

        MolWtOrg = [132.12;118.09;342.3];
        
        for ii = 1:N
            
            RH(1,ii,kkk,kkkk) = RHrange(ii);
            
            lognSD = 2;
            Ddry = DefineSizeDistribution(n,lognSD,k,mu,sigma); % Diameter of dry particle (um)
            Ddry = Ddry.*1e-6; % Diameter of dry particle (m)
            Ddry = sortrows(Ddry,1,'ascend'); % Sort dry diameters for plotting purpose
            Vdry = ((3.14/6)*(Ddry.^3)).*1000000000000; % Volume of dry particle (nl)
            
            Ddry = Ddry.*1000000; % Diameter of dry particle in um
            
            if isequal(core,'y')
                if BC_core == 1
                    DiaBC = 0.2+0.*Ddry; %[0.1;0.2;0.3;0.4]; % constant size of black carbon in um
                    Ddry = (Ddry.^3 + DiaBC.^3).^(1/3); % calculated by sum of the volumes
                elseif BC_core == 2
                    DiaBC = Ddry./1.1; % varying size of black carbon in um with const Dshell/Dcore
                    Ddry = (Ddry.^3 + DiaBC.^3).^(1/3); % calculated by sum of the volumes
                elseif BC_core == 3
                    shell_thickness = 0.5; % in um
                    DiaBC = Ddry-2*shell_thickness; % varying size of black carbon in um with const shell thickness (0.05 to 2.5 um)
                    Vdry = (3.14/6)*(Ddry.^3 - DiaBC.^3).*1e-6; % Volume of dry particle (nl)
                end
%                 Ddry = (Ddry.^3 + DiaBC.^3).^(1/3); % calculated by sum of the volumes
                ratio_Dshell_Dcore = Ddry./DiaBC;
                shell_thickness = (Ddry-DiaBC)./2;
            else
                DiaBC = 0;
            end
            
            massBC(:,1) = (((3.14/6)*(DiaBC.^3)).*(1e-12)).*1.8; % Mass of BC (g)

            pd_normalDry = fitdist(log(Ddry),'Normal'); % returns mu and sigma by fitting probability distribution to data
            
            SDdry = abs((2.302./((sqrt(2*3.14)).*((pd_normalDry.sigma)))).*(exp(-(((log(Ddry)-(pd_normalDry.mu)).^2)./(2.*((pd_normalDry.sigma)).^2)))));
            
            DiaMeanDry = exp(pd_normalDry.mu);
            DiaSigmaDry = exp(pd_normalDry.sigma);
            
            Solute1(:,ii) = ["NH42SO4";"NaCl"];
            
            Solute2(:,ii) = ["GlutaricAcid";"SuccinicAcid";"sucrose"];
            
            Solute(ii,:) = [Solute1(salt,ii),Solute2(org,ii)];
            
            massAS(:,ii,kkk,kkkk) = Vdry./(1/rhoSalt(salt) + OIR(:,:,kkk,kkkk)/rhoOrg(org)); % Mass of ammonium sulfate dry solute(ug)
            massOrg(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk).*OIR(:,:,kkk,kkkk); % Mass of organic dry solute(ug)
            
            vfAS(:,ii,kkk,kkkk) = (1/rhoSalt(salt))./((1/rhoSalt(salt)) + (OIR(:,:,kkk,kkkk)/rhoOrg(org)));
%             vfAS(:,ii,kkk,kkkk) = (massAS(:,ii,kkk,kkkk)/rhoSalt(salt))./Vdry;
            vfOrg(:,ii,kkk,kkkk) = (OIR(:,:,kkk,kkkk)/rhoOrg(org))./((1/rhoSalt(salt)) + (OIR(:,:,kkk,kkkk)/rhoOrg(org)));
            vfTot(:,ii,kkk,kkkk) = vfAS(:,ii,kkk,kkkk) + vfOrg(:,ii,kkk,kkkk);
            
            massSulfateAdsIso(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)./1.375; % Mass of sulfate from Adsorption Isotherm model(ug)
            
            massTotSolutes(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk).*(1+OIR(:,:,kkk,kkkk)); % Total mass of the solutes (ug)
            
            molarRatio(:,ii,kkk,kkkk) = (massAS(:,ii,kkk,kkkk)./MolWtSalt(salt))./(massOrg(:,ii,kkk,kkkk)./MolWtOrg(org)); % Molar ratio of salt to organic for use as an input in Adsorption Isotherm model (soluteratio in MolalityMixture)
            
            % %             aw = RH(1,ii,kkk,kkkk)/100; % for large particles
            
%             kappaSalt = 0.65; % AS
%             kappaSalt = 1.28; % NaCl
            kappaSalt = [0.65;1.28];
            
            kappaOrg = 0.1;
            
            kappa_vf(kkk) = vfAS(:,ii,kkk,kkkk).*kappaSalt(salt) + vfOrg(:,ii,kkk,kkkk).*kappaOrg;
            xt=size(Ddry);
            for pp = 1:xt
                
                [aw(pp,ii,kkk,kkkk),gf_kappa(pp,ii,kkk,kkkk),kappa_iter(pp,ii,kkk,kkkk)] = WaterActivity(kappa_vf(kkk),Ddry(pp),RH(1,ii,kkk,kkkk),STw,Vw,R,T); % kappa_iter should come out to be same as kappa_vf
%                 [aw(pp,ii,kkk,kkkk),STsol(pp,ii,kkk,kkkk),gf_kappa(pp,ii,kkk,kkkk),RH_iter(pp,ii,kkk,kkkk),kappa_iter(pp,ii,kkk,kkkk)] = WaterActivitySTsol(kappa_vf(kkk),Ddry(pp),RH(1,ii,kkk,kkkk),Vw,R,T); % valid for only binary AS solution

                [m(pp,:,ii,kkk,kkkk), Ix(pp,:,ii,kkk,kkkk), a(pp,:,ii,kkk,kkkk), f_Inf(pp,:,ii,kkk,kkkk), f_FS(pp,:,ii,kkk,kkkk), gamma_Inf(pp,:,ii,kkk,kkkk)] = MolalityMixture(aw(pp,ii,kkk,kkkk),molarRatio(1,1,kkk,kkkk),Solute(ii,:),2); % calculate molality, activity and coefficients of each component from Adsorption Isotherm model
                mTot(pp,ii,kkk,kkkk) = m(pp,1,ii,kkk,kkkk) + m(pp,2,ii,kkk,kkkk);
                
                [m1(pp,:,ii,kkk,kkkk), Ix1(pp,:,ii,kkk,kkkk), a1(pp,:,ii,kkk,kkkk), f_Inf1(pp,:,ii,kkk,kkkk), f_FS1(pp,:,ii,kkk,kkkk), gamma_Inf1(pp,:,ii,kkk,kkkk)] = MolalityMixture(aw(pp,ii,kkk,kkkk),1e+100,Solute(ii,:),2); % calculate molality, activity and coefficients of only the salt assuming no organic water uptake
                
            end
            
            massWaterAdsIso(:,ii,kkk,kkkk) = (1000.*(massAS(:,ii,kkk,kkkk)./MolWtSalt(salt) + massOrg(:,ii,kkk,kkkk)./MolWtOrg(org)))./mTot(:,ii,kkk,kkkk);
            
            massSol(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)+massOrg(:,ii,kkk,kkkk)+massWaterAdsIso(:,ii,kkk,kkkk);
            
            mfAS(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)./massSol(:,ii,kkk,kkkk); % Mass fraction of ammonium sulfate
            mfOrg(:,ii,kkk,kkkk) = massOrg(:,ii,kkk,kkkk)./massSol(:,ii,kkk,kkkk); % Mass fraction of organic
            mfWater(:,ii,kkk,kkkk) = massWaterAdsIso(:,ii,kkk,kkkk)./massSol(:,ii,kkk,kkkk);
            
            mfTot(:,ii,kkk,kkkk) = mfAS(:,ii,kkk,kkkk)+mfOrg(:,ii,kkk,kkkk)+mfWater(:,ii,kkk,kkkk); % check to verify it is 1
            
            for pp = 1:xt
                rhoSol(pp,ii,kkk,kkkk) = DensityMixture(mfAS(pp,ii,kkk,kkkk),mfOrg(pp,ii,kkk,kkkk),mfWater(pp,ii,kkk,kkkk),rhoSalt(salt),rhoOrg(org),ii); % calculate density of the solution from Laliberte model
            end
            
            VsolAdsIso(:,ii,kkk,kkkk) = massSol(:,ii,kkk,kkkk)./rhoSol(:,ii,kkk,kkkk); % Volume of solution in nl
            
            DwetAdsIso(:,ii,kkk,kkkk) = ((((VsolAdsIso(:,ii,kkk,kkkk)./1000000000000) + ((3.14/6)*((DiaBC*1e-6).^3)))./(3.14/6)).^(1/3)).*1000000; % Wet diameter in um
%             DwetAdsIso(:,ii,kkk,kkkk) = 0.3+0.*(((((VsolAdsIso(:,ii,kkk,kkkk)./1000000000000) + ((3.14/6)*((DiaBC*1e-6).^3)))./(3.14/6)).^(1/3)).*1000000); % const size for lensing effect
%             shell_thickness = 0.5;
%             DwetAdsIso(:,ii,kkk,kkkk) = shell_thickness.*2 + DiaBC;
%             shell_thickness(:,ii,kkk,kkkk) = (DwetAdsIso(:,ii,kkk,kkkk)-DiaBC)./2;
            ratio_Dshell_Dcore(:,ii,kkk,kkkk) = DwetAdsIso(:,ii,kkk,kkkk)./DiaBC;
            
            gf(:,ii,kkk,kkkk) = DwetAdsIso(:,ii,kkk,kkkk)./Ddry;
            kappa(:,ii,kkk,kkkk) = (gf(:,ii,kkk,kkkk).^3 - 1).*((1./aw(:,ii,kkk,kkkk)) - 1);
            
            pd_normal(:,ii,kkk,kkkk) = fitdist(log(DwetAdsIso(:,ii,kkk,kkkk)),'Normal'); % returns mu and sigma by fitting probability distribution to data
            DiaMean(:,ii,kkk,kkkk) = exp(pd_normal(:,ii,kkk,kkkk).mu);
            DiaSigma(:,ii,kkk,kkkk) = exp(pd_normal(:,ii,kkk,kkkk).sigma);
            
            %             SD(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*(log(pd_normal(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso(:,ii,kkk,kkkk))-log(pd_normal(:,ii,kkk,kkkk).mu)).^2)./(2.*(log(pd_normal(:,ii,kkk,kkkk).sigma)).^2)))));
            SD(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*((pd_normal(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso(:,ii,kkk,kkkk))-(pd_normal(:,ii,kkk,kkkk).mu)).^2)./(2.*((pd_normal(:,ii,kkk,kkkk).sigma)).^2)))));
            %%%% recalculation considering no water uptake by organics, yet just the mass of organics %%%%
            
            massWaterAdsIso1(:,ii,kkk,kkkk) = (1000.*massAS(:,ii,kkk,kkkk))./(MolWtSalt(salt).*m1(:,1,ii,kkk,kkkk)); % water uptake only by the salt
            massSol1(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)+massOrg(:,ii,kkk,kkkk)+massWaterAdsIso1(:,ii,kkk,kkkk);
            
            mfAS1(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)./massSol1(:,ii,kkk,kkkk); % Mass fraction of ammonium sulfate
            mfOrg1(:,ii,kkk,kkkk) = massOrg(:,ii,kkk,kkkk)./massSol1(:,ii,kkk,kkkk); % Mass fraction of organic
            mfWater1(:,ii,kkk,kkkk) = massWaterAdsIso1(:,ii,kkk,kkkk)./massSol1(:,ii,kkk,kkkk);
            
            mfTot1(:,ii,kkk,kkkk) = mfAS1(:,ii,kkk,kkkk)+mfOrg1(:,ii,kkk,kkkk)+mfWater1(:,ii,kkk,kkkk); % check to verify it is 1
            
            for pp = 1:xt
                rhoSol1(pp,ii,kkk,kkkk) = DensityMixture(mfAS1(pp,ii,kkk,kkkk),mfOrg1(pp,ii,kkk,kkkk),mfWater1(pp,ii,kkk,kkkk),rhoSalt(salt),rhoOrg(org),ii); % calculate density of the solution from Laliberte model
            end
            
            VsolAdsIso1(:,ii,kkk,kkkk) = massSol1(:,ii,kkk,kkkk)./rhoSol1(:,ii,kkk,kkkk); % Volume of solution in nl using actual density of mixture assuming organic contributes in density change
            
            DwetAdsIso1(:,ii,kkk,kkkk) = ((((VsolAdsIso1(:,ii,kkk,kkkk)./1000000000000) + ((3.14/6)*((DiaBC*1e-6).^3)))./(3.14/6)).^(1/3)).*1000000; % Wet diameter in um
%             DwetAdsIso1(:,ii,kkk,kkkk) = 0.3+0.*(((((VsolAdsIso1(:,ii,kkk,kkkk)./1000000000000) + ((3.14/6)*((DiaBC*1e-6).^3)))./(3.14/6)).^(1/3)).*1000000); % const size for lensing effect
%             DwetAdsIso1(:,ii,kkk,kkkk) = shell_thickness.*2 + DiaBC;

            gf1(:,ii,kkk,kkkk) = DwetAdsIso1(:,ii,kkk,kkkk)./Ddry;
            kappa1(:,ii,kkk,kkkk) = (gf1(:,ii,kkk,kkkk).^3 - 1).*((1./aw(:,ii,kkk,kkkk)) - 1);
            
            pd_normal1(:,ii,kkk,kkkk) = fitdist(log(DwetAdsIso1(:,ii,kkk,kkkk)),'Normal'); % returns mu and sigma by fitting probability distribution to data
            DiaMean1(:,ii,kkk,kkkk) = exp(pd_normal1(:,ii,kkk,kkkk).mu);
            DiaSigma1(:,ii,kkk,kkkk) = exp(pd_normal1(:,ii,kkk,kkkk).sigma);
            %             SD1(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*(log(pd_normal1(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso1(:,ii,kkk,kkkk))-log(pd_normal1(:,ii,kkk,kkkk).mu)).^2)./(2.*(log(pd_normal1(:,ii,kkk,kkkk).sigma)).^2)))));
            SD1(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*((pd_normal1(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso1(:,ii,kkk,kkkk))-(pd_normal1(:,ii,kkk,kkkk).mu)).^2)./(2.*((pd_normal1(:,ii,kkk,kkkk).sigma)).^2)))));
%                         DiaMeanError(:,ii,kkk,kkkk) = ((DiaMean(:,ii,kkk,kkkk) - DiaMean1(:,ii,kkk,kkkk))./DiaMean(:,ii,kkk,kkkk)).*100; % error(%) between the two water uptake calculations
            
            % To find Dwet_kappa, use eqn 2 from [PK07] --> 1/aw = 1 + kappa(Vdry/Vwater)
            Vwater_kappa(:,ii,kkk,kkkk) = (kappa_iter(:,ii,kkk,kkkk).*Vdry)./((1./aw(:,ii,kkk,kkkk)) - 1);
            Vsol_kappa(:,ii,kkk,kkkk) = Vwater_kappa(:,ii,kkk,kkkk) + Vdry;
            Dwet_kappa(:,ii,kkk,kkkk) = ((((Vsol_kappa(:,ii,kkk,kkkk)./1000000000000) + ((3.14/6)*((DiaBC*1e-6).^3)))./(3.14/6)).^(1/3)).*1000000; % Wet diameter in um
            
            pd_normal2(:,ii,kkk,kkkk) = fitdist(log(Dwet_kappa(:,ii,kkk,kkkk)),'Normal'); % returns mu and sigma by fitting probability distribution to data
            DiaMean2(:,ii,kkk,kkkk) = exp(pd_normal2(:,ii,kkk,kkkk).mu);
            DiaSigma2(:,ii,kkk,kkkk) = exp(pd_normal2(:,ii,kkk,kkkk).sigma);
            SD2(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*((pd_normal2(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(Dwet_kappa(:,ii,kkk,kkkk))-(pd_normal2(:,ii,kkk,kkkk).mu)).^2)./(2.*((pd_normal2(:,ii,kkk,kkkk).sigma)).^2)))));
            
            xSizeDry = (pi.*Ddry)./(lambda(1,1,1,kkkk)*0.001); % size parameter
            xSize(:,ii,kkk,kkkk) = (pi.*DwetAdsIso(:,ii,kkk,kkkk))./(lambda(1,1,1,kkkk)*0.001); % size parameter
            xSize1(:,ii,kkk,kkkk) = (pi.*DwetAdsIso1(:,ii,kkk,kkkk))./(lambda(1,1,1,kkkk)*0.001); % size parameter
            xSize_kappa(:,ii,kkk,kkkk) = (pi.*Dwet_kappa(:,ii,kkk,kkkk))./(lambda(1,1,1,kkkk)*0.001); % size parameter
            
            massWater_kappa(:,ii,kkk,kkkk) = Vwater_kappa(:,ii,kkk,kkkk).*0.998;
            massSol2(:,ii,kkk,kkkk) = massTotSolutes(:,ii,kkk,kkkk) + massWater_kappa(:,ii,kkk,kkkk);
            rhoSol2(:,ii,kkk,kkkk) = massSol2(:,ii,kkk,kkkk)./Vsol_kappa(:,ii,kkk,kkkk);
            
            mfAS2(:,ii,kkk,kkkk) = massAS(:,ii,kkk,kkkk)./massSol2(:,ii,kkk,kkkk); % Mass fraction of ammonium sulfate
            mfOrg2(:,ii,kkk,kkkk) = massOrg(:,ii,kkk,kkkk)./massSol2(:,ii,kkk,kkkk); % Mass fraction of organic
            mfWater2(:,ii,kkk,kkkk) = massWater_kappa(:,ii,kkk,kkkk)./massSol2(:,ii,kkk,kkkk);
            
            mfTot2(:,ii,kkk,kkkk) = mfAS2(:,ii,kkk,kkkk)+mfOrg2(:,ii,kkk,kkkk)+mfWater2(:,ii,kkk,kkkk); % check to verify it is 1
            
            for pp = 1:xt
                rhoSol2_recheck(pp,ii,kkk,kkkk) = DensityMixture(mfAS2(pp,ii,kkk,kkkk),mfOrg2(pp,ii,kkk,kkkk),mfWater2(pp,ii,kkk,kkkk),rhoSalt(salt),rhoOrg(org),ii); % calculate density of the solution from Laliberte model to compare with rhoSol2
            end
            
%             RIsaltR = 1.5298; % refractive index of the dry salt(AS) at all visible wavelengths (only the real part)
            %             RIsaltIrange = [0.005,0.004,0.005,0.013,0.012,0.018,0.023,0.027,0.0339];
%             RIsaltR = 1.553; % refractive index of NaCl at 550 nm
            RIsaltR = [1.5298;1.553];
            %             RIsaltI(1,1,1,kkkk) = RIsaltIrange(kkkk); % refractive index of the dry salt(AS) (only the imaginary part)
            
            RIorgR = [1.42;1.45;1.5]; % refractive index of the dry organic at all visible wavelengths (only the real part)
            RIorg = RIorgR(org) + 0.0i;
            
            RIwaterRrange = 1.333; % [1.343,1.339,1.337,1.335,1.333,1.332,1.331,1.331,1.330]; % refractive index of pure water at wavelengths 350nm:50:750nm (only the real part)
            RIwaterR(1,1,1,kkkk) = RIwaterRrange(kkkk);
            
            RIreal(:,ii,kkk,kkkk) = rhoSol(:,ii,kkk,kkkk).*((mfAS(:,ii,kkk,kkkk)./rhoSalt(salt)).*RIsaltR(salt) + (mfOrg(:,ii,kkk,kkkk)./rhoOrg(org)).*RIorgR(org) + (mfWater(:,ii,kkk,kkkk)./0.998).*RIwaterR(1,1,1,kkkk)); % refractive index of the solution (only the real part)
            RIreal1(:,ii,kkk,kkkk) = rhoSol1(:,ii,kkk,kkkk).*((mfAS1(:,ii,kkk,kkkk)./rhoSalt(salt)).*RIsaltR(salt) + (mfOrg1(:,ii,kkk,kkkk)./rhoOrg(org)).*RIorgR(org) + (mfWater1(:,ii,kkk,kkkk)./0.998).*RIwaterR(1,1,1,kkkk)); % refractive index of the solution (only the real part)
            RIreal2(:,ii,kkk,kkkk) = rhoSol2(:,ii,kkk,kkkk).*((mfAS2(:,ii,kkk,kkkk)./rhoSalt(salt)).*RIsaltR(salt) + (mfOrg2(:,ii,kkk,kkkk)./rhoOrg(org)).*RIorgR(org) + (mfWater2(:,ii,kkk,kkkk)./0.998).*RIwaterR(1,1,1,kkkk)); % refractive index of the solution (only the real part)
            
            RIerror1(:,ii,kkk,kkkk) = ((RIreal(:,ii,kkk,kkkk) - RIreal1(:,ii,kkk,kkkk))./RIreal(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations            RIsalt = RIsaltR + 0.0i;
            RIerror2(:,ii,kkk,kkkk) = ((RIreal(:,ii,kkk,kkkk) - RIreal2(:,ii,kkk,kkkk))./RIreal(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations            RIsalt = RIsaltR + 0.0i;
            
%             RI(:,ii,kkk,kkkk) = 1.38+0.*RIreal(:,ii,kkk,kkkk) + 0.0i; % const RI to check lensing effect
%             RI1(:,ii,kkk,kkkk) = 1.38+0.*RIreal(:,ii,kkk,kkkk) + 0.0i; % const RI to check lensing effect
            RI(:,ii,kkk,kkkk) = RIreal(:,ii,kkk,kkkk) + 0.0i; % *RIsaltI(1,1,1,kkkk);
            RI1(:,ii,kkk,kkkk) = RIreal1(:,ii,kkk,kkkk) + 0.0i; % *RIsaltI(1,1,1,kkkk);
            RI2(:,ii,kkk,kkkk) = RIreal2(:,ii,kkk,kkkk) + 0.0i; % *RIsaltI(1,1,1,kkkk);
            
            if isequal(core,'y')
                %                 DiaBC = 0.08 + 0*Ddry./2; % size of black carbon in um
                %                 Ddry = (Ddry.^3 + DiaBC.^3).^(1/3); % calculated by sum of the volumes, not diameter
                ySize(:,ii,kkk,kkkk) = (pi.*DiaBC)./(lambda(1,1,1,kkkk)*0.001); % size parameter of black carbon
                %                                 RIbc(1,ii,kkk,kkkk) = 1.82 + 0.74i; % RI of black carbon at 550 nm
                RIbc(1,ii,kkk,kkkk) = 1.95 + 0.79i; % RI of black carbon at 550 nm (Ref: Bond and Bergstrom, 2006))
                %                 RIbcRange = 1.75+0.44i; % [1.75+0.465i,1.75+0.46i,1.75+0.455i,1.75+0.45i,1.75+0.44i,1.75+0.435i,1.75+0.435i,1.75+0.43i,1.75+0.43i];
                %                 RIbc(1,ii,kkk,kkkk) = RIbcRange(kkkk); % RI of black carbon
                
                for pp = 1:xt
                    
                    [qscaBC(pp,ii,kkk,kkkk),qabsBC(pp,ii,kkk,kkkk)] = Mie(ySize(pp,ii,kkk,kkkk),RIbc(:,ii,kkk,kkkk)); % Mie scattering subroutine for the bare core (no coating)
                    [qext(pp,ii,kkk,kkkk),qsca(pp,ii,kkk,kkkk),qabs(pp,ii,kkk,kkkk),qback(pp,ii,kkk,kkkk),gsca(pp,ii,kkk,kkkk),bratio(pp,ii,kkk,kkkk)]=Miecoated(xSize(pp,ii,kkk,kkkk),ySize(pp,ii,kkk,kkkk),RI(pp,ii,kkk,kkkk),RIbc(:,ii,kkk,kkkk),opt); % Mie scattering subroutine for coated droplets
                    [qext1(pp,ii,kkk,kkkk),qsca1(pp,ii,kkk,kkkk),qabs1(pp,ii,kkk,kkkk),qback1(pp,ii,kkk,kkkk),gsca1(pp,ii,kkk,kkkk),bratio1(pp,ii,kkk,kkkk)]=Miecoated(xSize1(pp,ii,kkk,kkkk),ySize(pp,ii,kkk,kkkk),RI1(pp,ii,kkk,kkkk),RIbc(:,ii,kkk,kkkk),opt); % Mie scattering subroutine for coated droplets
                    [qext2(pp,ii,kkk,kkkk),qsca2(pp,ii,kkk,kkkk),qabs2(pp,ii,kkk,kkkk),qback2(pp,ii,kkk,kkkk),gsca2(pp,ii,kkk,kkkk),bratio2(pp,ii,kkk,kkkk)]=Miecoated(xSize_kappa(pp,ii,kkk,kkkk),ySize(pp,ii,kkk,kkkk),RI2(pp,ii,kkk,kkkk),RIbc(:,ii,kkk,kkkk),opt); % Mie scattering subroutine for coated droplets
                    ext(pp,ii,kkk,kkkk) = qext(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);
                    ext1(pp,ii,kkk,kkkk) = qext1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    ext2(pp,ii,kkk,kkkk) = qext2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    backscat(pp,ii,kkk,kkkk) = qback(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);
                    backscat1(pp,ii,kkk,kkkk) = qback1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    backscat2(pp,ii,kkk,kkkk) = qback2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    
                    scatBC(pp,ii,kkk,kkkk) = qscaBC(pp,ii,kkk,kkkk).*(pi/4).*(DiaBC(pp).^2);
                    absorpBC(pp,ii,kkk,kkkk) = qabsBC(pp,ii,kkk,kkkk).*(pi/4).*(DiaBC(pp).^2);
                    
                    scat(pp,ii,kkk,kkkk) = qsca(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);  % Scattering cross-section of each particle size (scattering efficiency*surface area of each particle diameter) (um^2)
                    scat1(pp,ii,kkk,kkkk) = qsca1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    scat2(pp,ii,kkk,kkkk) = qsca2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    
                    scatEnhance(pp,ii,kkk,kkkk) = scat(pp,ii,kkk,kkkk)./scatBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    scatEnhance1(pp,ii,kkk,kkkk) = scat1(pp,ii,kkk,kkkk)./scatBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    scatEnhance2(pp,ii,kkk,kkkk) = scat2(pp,ii,kkk,kkkk)./scatBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    
                    absorp(pp,ii,kkk,kkkk) = qabs(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);
                    absorp1(pp,ii,kkk,kkkk) = qabs1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    absorp2(pp,ii,kkk,kkkk) = qabs2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    
                    MACbc(pp,ii,kkk,kkkk) = (absorp(pp,ii,kkk,kkkk).*(1e-12))./massBC(pp); % normalized mass abs c/s in m2/g, also termed as MAC in some papers
                    MACcore(pp,ii,kkk,kkkk) = (absorpBC(pp,ii,kkk,kkkk).*(1e-12))./massBC(pp); % normalized mass abs c/s in m2/g
                                        
                    absorpEnhance(pp,ii,kkk,kkkk) = absorp(pp,ii,kkk,kkkk)./absorpBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    absorpEnhance1(pp,ii,kkk,kkkk) = absorp1(pp,ii,kkk,kkkk)./absorpBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    absorpEnhance2(pp,ii,kkk,kkkk) = absorp2(pp,ii,kkk,kkkk)./absorpBC(pp,ii,kkk,kkkk); % aerosol c/s absorption enhancement
                    
                    particles = n;
                    %                     particles = 1e+4; % #concentration (#/m^3 volume of air)
                    
                    nsteps = 2*angleNo - 1;
                    [angularScatPF(:,:,pp,ii,kkk,kkkk),s1(pp,:,ii,kkk,kkkk),s2(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI(pp,ii,kkk,kkkk), xSize(pp), nsteps);
                    [angularScatPF1(:,:,pp,ii,kkk,kkkk),s1_1(pp,:,ii,kkk,kkkk),s2_1(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI1(pp,ii,kkk,kkkk), xSize1(pp), nsteps);
                    [angularScatPF2(:,:,pp,ii,kkk,kkkk),s1_2(pp,:,ii,kkk,kkkk),s2_2(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI2(pp,ii,kkk,kkkk), xSize_kappa(pp), nsteps);
                    ScatPF(pp,:,ii,kkk,kkkk) = ((DwetAdsIso(pp,ii,kkk,kkkk).^2).*(abs(s1(pp,:,ii,kkk,kkkk)).^2+abs(s2(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca(pp,ii,kkk,kkkk));
                    ScatPF1(pp,:,ii,kkk,kkkk) = ((DwetAdsIso1(pp,ii,kkk,kkkk).^2).*(abs(s1_1(pp,:,ii,kkk,kkkk)).^2+abs(s2_1(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca1(pp,ii,kkk,kkkk));
                    ScatPF2(pp,:,ii,kkk,kkkk) = ((Dwet_kappa(pp,ii,kkk,kkkk).^2).*(abs(s1_2(pp,:,ii,kkk,kkkk)).^2+abs(s2_2(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca2(pp,ii,kkk,kkkk));
                    
                end
                
                TotExt(:,ii,kkk,kkkk) = sum(ext(:,ii,kkk,kkkk))*1e-9;
                TotExt1(:,ii,kkk,kkkk) = sum(ext1(:,ii,kkk,kkkk))*1e-9;
                TotExt2(:,ii,kkk,kkkk) = sum(ext2(:,ii,kkk,kkkk))*1e-9;
                
                TotScatBC(:,ii,kkk,kkkk) = sum(scatBC(:,ii,kkk,kkkk))*1e-9;
                TotScat(:,ii,kkk,kkkk) = sum(scat(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotScat1(:,ii,kkk,kkkk) = sum(scat1(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotScat2(:,ii,kkk,kkkk) = sum(scat2(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                
                ScatCalcError1(:,ii,kkk,kkkk) = ((TotScat(:,ii,kkk,kkkk) - TotScat1(:,ii,kkk,kkkk))./TotScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                ScatCalcError2(:,ii,kkk,kkkk) = ((TotScat(:,ii,kkk,kkkk) - TotScat2(:,ii,kkk,kkkk))./TotScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
                TotScatEnhance(:,ii,kkk,kkkk) = TotScat(:,ii,kkk,kkkk)./TotScatBC(:,ii,kkk,kkkk); % aerosol population absorption enhancement
                TotScatEnhance1(:,ii,kkk,kkkk) = TotScat1(:,ii,kkk,kkkk)./TotScatBC(:,ii,kkk,kkkk);
                TotScatEnhance2(:,ii,kkk,kkkk) = TotScat2(:,ii,kkk,kkkk)./TotScatBC(:,ii,kkk,kkkk);
                
                TotBackScat(:,ii,kkk,kkkk) = sum(backscat(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotBackScat1(:,ii,kkk,kkkk) = sum(backscat1(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotBackScat2(:,ii,kkk,kkkk) = sum(backscat2(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                
                BackScatCalcError1(:,ii,kkk,kkkk) = ((TotBackScat(:,ii,kkk,kkkk) - TotBackScat1(:,ii,kkk,kkkk))./TotBackScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                BackScatCalcError2(:,ii,kkk,kkkk) = ((TotBackScat(:,ii,kkk,kkkk) - TotBackScat2(:,ii,kkk,kkkk))./TotBackScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
                TotAbsBC(:,ii,kkk,kkkk) = sum(absorpBC(:,ii,kkk,kkkk))*1e-9;
                TotAbs(:,ii,kkk,kkkk)=sum(absorp(:,ii,kkk,kkkk))*1e-9;
%                 MACbc(:,ii,kkk,kkkk) = (TotAbs(:,ii,kkk,kkkk).*(1e-12))./massBC(1); % normalized mass abs coeff of BC alone in m2/g
                TotAbs1(:,ii,kkk,kkkk)=sum(absorp1(:,ii,kkk,kkkk))*1e-9;
                TotAbs2(:,ii,kkk,kkkk)=sum(absorp2(:,ii,kkk,kkkk))*1e-9;
                
                MACtotal(:,ii,kkk,kkkk) = (TotAbs(:,ii,kkk,kkkk).*(1e-9))./(sum(massBC)); % normalized mass abs c/s in m2/g

%                MACtotal(:,ii,kkk,kkkk) = (sum(absorp(:,ii,kkk,kkkk)).*(1e-12))./(sum(massBC(:))); % normalized mass abs c/s in m2/g

                
                AbsCalcError1(:,ii,kkk,kkkk) = ((TotAbs(:,ii,kkk,kkkk) - TotAbs1(:,ii,kkk,kkkk))./TotAbs(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                AbsCalcError2(:,ii,kkk,kkkk) = ((TotAbs(:,ii,kkk,kkkk) - TotAbs2(:,ii,kkk,kkkk))./TotAbs(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
                TotAbsEnhance(:,ii,kkk,kkkk) = TotAbs(:,ii,kkk,kkkk)./TotAbsBC(:,ii,kkk,kkkk); % aerosol population absorption enhancement
                TotAbsEnhance1(:,ii,kkk,kkkk) = TotAbs1(:,ii,kkk,kkkk)./TotAbsBC(:,ii,kkk,kkkk);
                TotAbsEnhance2(:,ii,kkk,kkkk) = TotAbs2(:,ii,kkk,kkkk)./TotAbsBC(:,ii,kkk,kkkk);
                
                SSAbc(:,ii,kkk,kkkk) = TotScatBC(:,ii,kkk,kkkk)./(TotScatBC(:,ii,kkk,kkkk)+TotAbsBC(:,ii,kkk,kkkk)); % single scattering albedo
                SSA(:,ii,kkk,kkkk) = TotScat(:,ii,kkk,kkkk)./TotExt(:,ii,kkk,kkkk); % single scattering albedo
                SSA1(:,ii,kkk,kkkk) = TotScat1(:,ii,kkk,kkkk)./TotExt1(:,ii,kkk,kkkk); % single scattering albedo
                SSA2(:,ii,kkk,kkkk) = TotScat2(:,ii,kkk,kkkk)./TotExt2(:,ii,kkk,kkkk); % single scattering albedo
                
                SSAerror1(:,ii,kkk,kkkk) = ((SSA(:,ii,kkk,kkkk) - SSA1(:,ii,kkk,kkkk))./SSA(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                SSAerror2(:,ii,kkk,kkkk) = ((SSA(:,ii,kkk,kkkk) - SSA2(:,ii,kkk,kkkk))./SSA(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
            else
                
                for pp = 1:xt
                    [qext(pp,ii,kkk,kkkk),qsca(pp,ii,kkk,kkkk),qback(pp,ii,kkk,kkkk),gsca(pp,ii,kkk,kkkk),bratio(pp,ii,kkk,kkkk)]=Mie(xSize(pp,ii,kkk,kkkk),RI(pp,ii,kkk,kkkk)); % Mie scattering subroutine for homogeneous organic-inorganic droplets (unitless)
                    
                    %                     [s1_bh(pp,:,ii,kkk,kkkk),s2_bh(pp,:,ii,kkk,kkkk),gsca_bh(pp,ii,kkk,kkkk)]=bhmie(xSize(pp,ii,kkk,kkkk),RI(:,ii,kkk,kkkk),angleNo);
                    %                     p_s1s2_bh(pp,:,ii,kkk,kkkk) = ((DwetAdsIso(pp,ii,kkk,kkkk).^2).*(abs(s1_bh(pp,:,ii,kkk,kkkk)).^2+abs(s2_bh(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca(pp,ii,kkk,kkkk));
                    %                     p_g(pp,:,ii,kkk,kkkk) = (1/(4*3.14)).*((1-(gsca(pp,ii,kkk,kkkk).^2))./(1+(gsca(pp,ii,kkk,kkkk).^2)-(2.*gsca(pp,ii,kkk,kkkk).*cos(theta*(3.14/180)))));
                    %                     bratio_g(pp,ii,kkk,kkkk) = ((1-gsca(pp,ii,kkk,kkkk))./(2.*gsca(pp,ii,kkk,kkkk))).*(((1-gsca(pp,ii,kkk,kkkk))./sqrt(1+(gsca(pp,ii,kkk,kkkk).^2)))-1);
                    
                    [qext1(pp,ii,kkk,kkkk),qsca1(pp,ii,kkk,kkkk),qback1(pp,ii,kkk,kkkk),gsca1(pp,ii,kkk,kkkk),bratio1(pp,ii,kkk,kkkk)]=Mie(xSize1(pp,ii,kkk,kkkk),RI1(pp,ii,kkk,kkkk)); % Mie scattering subroutine for same droplets but water uptake only by inorganic
                    [qext2(pp,ii,kkk,kkkk),qsca2(pp,ii,kkk,kkkk),qback2(pp,ii,kkk,kkkk),gsca2(pp,ii,kkk,kkkk),bratio2(pp,ii,kkk,kkkk)]=Mie(xSize_kappa(pp,ii,kkk,kkkk),RI2(pp,ii,kkk,kkkk)); % Mie scattering subroutine for same droplets but water uptake using hygroscopicity paramaters (kappa)
                    ext(pp,ii,kkk,kkkk) = qext(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);
                    ext1(pp,ii,kkk,kkkk) = qext1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    ext2(pp,ii,kkk,kkkk) = qext2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    scat(pp,ii,kkk,kkkk) = qsca(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);  % Scattering cross-section of each particle size (scattering efficiency*surface area of each particle diameter) (um^2)
                    scat1(pp,ii,kkk,kkkk) = qsca1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    scat2(pp,ii,kkk,kkkk) = qsca2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    backscat(pp,ii,kkk,kkkk) = qback(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso(pp,ii,kkk,kkkk).^2);
                    backscat1(pp,ii,kkk,kkkk) = qback1(pp,ii,kkk,kkkk).*(pi/4).*(DwetAdsIso1(pp,ii,kkk,kkkk).^2);
                    backscat2(pp,ii,kkk,kkkk) = qback2(pp,ii,kkk,kkkk).*(pi/4).*(Dwet_kappa(pp,ii,kkk,kkkk).^2);
                    
                    %                     backscat_g(pp,ii,kkk,kkkk) = scat(pp,ii,kkk,kkkk).*bratio_g(pp,ii,kkk,kkkk);
                    
                    particles = n;
                    %                     particles = 1e+4; % #concentration (#/m^3 volume of air)
                    
                    nsteps = 2*angleNo - 1;
                    [angularScatPF(:,:,pp,ii,kkk,kkkk),s1(pp,:,ii,kkk,kkkk),s2(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI(pp,ii,kkk,kkkk), xSize(pp), nsteps);
                    [angularScatPF1(:,:,pp,ii,kkk,kkkk),s1_1(pp,:,ii,kkk,kkkk),s2_1(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI1(pp,ii,kkk,kkkk), xSize1(pp), nsteps);
                    [angularScatPF2(:,:,pp,ii,kkk,kkkk),s1_2(pp,:,ii,kkk,kkkk),s2_2(pp,:,ii,kkk,kkkk)] = Mie_tetascan(RI2(pp,ii,kkk,kkkk), xSize_kappa(pp), nsteps);
                    ScatPF(pp,:,ii,kkk,kkkk) = ((DwetAdsIso(pp,ii,kkk,kkkk).^2).*(abs(s1(pp,:,ii,kkk,kkkk)).^2+abs(s2(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca(pp,ii,kkk,kkkk));
                    ScatPF1(pp,:,ii,kkk,kkkk) = ((DwetAdsIso1(pp,ii,kkk,kkkk).^2).*(abs(s1_1(pp,:,ii,kkk,kkkk)).^2+abs(s2_1(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca1(pp,ii,kkk,kkkk));
                    ScatPF2(pp,:,ii,kkk,kkkk) = ((Dwet_kappa(pp,ii,kkk,kkkk).^2).*(abs(s1_2(pp,:,ii,kkk,kkkk)).^2+abs(s2_2(pp,:,ii,kkk,kkkk)).^2))./(8.*qsca2(pp,ii,kkk,kkkk));
                    
                end
                
                TotExt(:,ii,kkk,kkkk) = sum(ext(:,ii,kkk,kkkk))*1e-9;
                TotExt1(:,ii,kkk,kkkk) = sum(ext1(:,ii,kkk,kkkk))*1e-9;
                TotExt2(:,ii,kkk,kkkk) = sum(ext2(:,ii,kkk,kkkk))*1e-9;
                
                TotScat(:,ii,kkk,kkkk) = sum(scat(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotScat1(:,ii,kkk,kkkk) = sum(scat1(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotScat2(:,ii,kkk,kkkk) = sum(scat2(:,ii,kkk,kkkk))*1e-9; % Total scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                
                ScatCalcError1(:,ii,kkk,kkkk) = ((TotScat(:,ii,kkk,kkkk) - TotScat1(:,ii,kkk,kkkk))./TotScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                ScatCalcError2(:,ii,kkk,kkkk) = ((TotScat(:,ii,kkk,kkkk) - TotScat2(:,ii,kkk,kkkk))./TotScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
                TotBackScat(:,ii,kkk,kkkk) = sum(backscat(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotBackScat1(:,ii,kkk,kkkk) = sum(backscat1(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                TotBackScat2(:,ii,kkk,kkkk) = sum(backscat2(:,ii,kkk,kkkk))*1e-9; % Total back-scattering coefficient (#concentration of each size bin*scattering cross-section of each particle size) (1e-18*um^2 / 1e-9*m^3 = km^-1)
                
                BackScatCalcError1(:,ii,kkk,kkkk) = ((TotBackScat(:,ii,kkk,kkkk) - TotBackScat1(:,ii,kkk,kkkk))./TotBackScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                BackScatCalcError2(:,ii,kkk,kkkk) = ((TotBackScat(:,ii,kkk,kkkk) - TotBackScat2(:,ii,kkk,kkkk))./TotBackScat(:,ii,kkk,kkkk)).*100; % error(%) between different water uptake calculations
                
                TotAbs(:,ii,kkk,kkkk)=TotExt(:,ii,kkk,kkkk)-TotScat(:,ii,kkk,kkkk);
                TotAbs1(:,ii,kkk,kkkk)=TotExt1(:,ii,kkk,kkkk)-TotScat1(:,ii,kkk,kkkk);
                TotAbs2(:,ii,kkk,kkkk)=TotExt2(:,ii,kkk,kkkk)-TotScat2(:,ii,kkk,kkkk);
                
                SSA(:,ii,kkk,kkkk) = TotScat(:,ii,kkk,kkkk)./TotExt(:,ii,kkk,kkkk); % single scattering albedo
                SSA1(:,ii,kkk,kkkk) = TotScat1(:,ii,kkk,kkkk)./TotExt1(:,ii,kkk,kkkk); % single scattering albedo
                SSA2(:,ii,kkk,kkkk) = TotScat2(:,ii,kkk,kkkk)./TotExt2(:,ii,kkk,kkkk); % single scattering albedo
                
            end
            
            RHid = sprintf(' RH = %d',RH(1,ii,kkk,kkkk));
            OIRid = sprintf(' OIR = %d',OIR(1,1,kkk,kkkk));
            lambdaID = sprintf('Wavelength = %d nm',lambda(1,1,1,kkkk));
            Solute1id = string(Solute1(salt,ii));
            Solute2id = string(Solute2(org,ii));
            BCid = sprintf(' DiaBC = %d',DiaBC(1));
            
        end
        
        % Normalized mean square error between the two water uptake calculations:
        % AdsIso model and No org water uptake model
        Norm_MSE_Dia1(:,:,kkk,kkkk) = (sum(((DiaMean(:,:,kkk,kkkk) - DiaMean1(:,:,kkk,kkkk))./DiaMean(:,:,kkk,kkkk)).^2)/N(1));
        
        Norm_MSE_RI1(:,:,kkk,kkkk) = (sum(((RIreal(:,:,kkk,kkkk) - RIreal1(:,:,kkk,kkkk))./RIreal(:,:,kkk,kkkk)).^2)/N(1));
        
        Norm_MSE_Scat1(:,:,kkk,kkkk) = (sum(((TotScat(:,:,kkk,kkkk) - TotScat1(:,:,kkk,kkkk))./TotScat(:,:,kkk,kkkk)).^2)/N(1));
        Norm_MSE_BackScat1(:,:,kkk,kkkk) = (sum(((TotBackScat(:,:,kkk,kkkk) - TotBackScat1(:,:,kkk,kkkk))./TotBackScat(:,:,kkk,kkkk)).^2)/N(1)); %
        Norm_MSE_Abs1(:,:,kkk,kkkk) = (sum(((TotAbs(:,:,kkk,kkkk) - TotAbs1(:,:,kkk,kkkk))./TotAbs(:,:,kkk,kkkk)).^2)/N(1));
        Norm_MSE_SSA1(:,:,kkk,kkkk) = (sum(((SSA(:,:,kkk,kkkk) - SSA1(:,:,kkk,kkkk))./SSA(:,:,kkk,kkkk)).^2)/N(1));
        
        % Normalized mean square error between the two water uptake calculations:
        % AdsIso model and Hygroscopicity kappa model
        Norm_MSE_Dia2(:,:,kkk,kkkk) = (sum(((DiaMean(:,:,kkk,kkkk) - DiaMean2(:,:,kkk,kkkk))./DiaMean(:,:,kkk,kkkk)).^2)/N(1));
        
        Norm_MSE_RI2(:,:,kkk,kkkk) = (sum(((RIreal(:,:,kkk,kkkk) - RIreal2(:,:,kkk,kkkk))./RIreal(:,:,kkk,kkkk)).^2)/N(1));
        
        Norm_MSE_Scat2(:,:,kkk,kkkk) = (sum(((TotScat(:,:,kkk,kkkk) - TotScat2(:,:,kkk,kkkk))./TotScat(:,:,kkk,kkkk)).^2)/N(1));
        Norm_MSE_BackScat2(:,:,kkk,kkkk) = (sum(((TotBackScat(:,:,kkk,kkkk) - TotBackScat2(:,:,kkk,kkkk))./TotBackScat(:,:,kkk,kkkk)).^2)/N(1)); %
        Norm_MSE_Abs2(:,:,kkk,kkkk) = (sum(((TotAbs(:,:,kkk,kkkk) - TotAbs1(:,:,kkk,kkkk))./TotAbs(:,:,kkk,kkkk)).^2)/N(1));
        Norm_MSE_SSA2(:,:,kkk,kkkk) = (sum(((SSA(:,:,kkk,kkkk) - SSA1(:,:,kkk,kkkk))./SSA(:,:,kkk,kkkk)).^2)/N(1));
        
        %%%%%%%% -- For each case, write all data in log files of individual -- %%%%%%%%
        %%%%%%%% -- wavelength, with individual OIR tab in each sheet -- %%%%%%%%
        %%%%%%%% -- having data of particle property/variable versus RH.  -- %%%%%%%%
        %%%%%%%% -- writematrix(A,filename) writes to a file with the -- %%%%%%%%
        %%%%%%%% -- name and extension specified by filename. -- %%%%%%%%
        %%%%%%%% -- But first, create a folder with date. -- %%%%%%%%
        %%%%%%%% -- This folder directory will store log files with -- %%%%%%%%
        %%%%%%%% -- each solutes' name and individual wavelength. -- %%%%%%%%
        
        foldername = sprintf('%s/',now);
        folderpath = strcat('./Water-Uptake_Salt-Org_Opt-properties/',foldername);
        warning("off")
        status = mkdir(folderpath);
        fprintf(folderpath,'\n');
        
    end
        
%         if isequal(core,'y')
%             ToptProp = table(RHrange,SSAbc,SSA,TotAbsBC,TotAbs,TotScatBC,TotScat,TotAbsEnhance,TotScatEnhance);
%             filename16 = strcat(folderpath,'/','OptProperties','-',Solute1id,'+', Solute2id,'-',BCid,'.xlsx');
%             writetable(ToptProp,filename16);
%         end
        
end

matfile = strcat(folderpath,'/',Solute1id,'+', Solute2id,'-',core,'.mat');
save(matfile);

runtimeseconds = toc(starttimer);

        % %             TabsEnhance = table(TotAbsEnhance,TotAbsEnhance1,TotAbsEnhance2);
        % %             filename15 = strcat(folderpath,'/','AbsEnhance','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %             writetable(TabsEnhance,filename15);
        
        % %         Taw = table(aw);
        % %         TSD = table(SDdry,DwetAdsIso1,SD1,DwetAdsIso,SD,Dwet_kappa,SD2);
        % %         TGF = table(gf1,gf,gf_kappa);
        % %         TKappa = table(kappa1,kappa,kappa_iter);
        % %         TScatError = table(ScatCalcError1,ScatCalcError2);
        % %         TAbsError = table(AbsCalcError1,AbsCalcError2);
        % %         TSSAerror = table(SSAerror1,SSAerror2);
        
        % %         filename7 = strcat(folderpath,'/','aw','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename8 = strcat(folderpath,'/','SD','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename9 = strcat(folderpath,'/','GF','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename11 = strcat(folderpath,'/','Kappa','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename12 = strcat(folderpath,'/','ScatError','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename13 = strcat(folderpath,'/','AbsError','-',Solute1id,'+', Solute2id,'-','.xlsx');
        % %         filename14 = strcat(folderpath,'/','SSAerror','-',Solute1id,'+', Solute2id,'-','.xlsx');
        
        % %         writetable(Taw,filename7);
        % %         writetable(TSD,filename8);
        % %         writetable(TGF,filename9);
        % %         writetable(TKappa,filename11);
        % %         writetable(TScatError,filename12);
        % %         writetable(TAbsError,filename13);
        % %         writetable(TSSAerror,filename14);
        
%     Ddry = [1e-9:8e-10:1e-8,1.01e-8:5e-9:1e-7,1.01e-7:1e-8:1e-6,1.01e-6:1e-8:1e-5,1.01e-5:1e-7:1e-4]'; % Diameter of dry particle (m) in ascending order

% %         lognSD=2;
% %         Ddry = DefineSizeDistribution(n,lognSD,k,mu,sigma); % Diameter of dry particle (um)
% %         Ddry = Ddry.*1e-6; % Diameter of dry particle (m)
% %         %         [1e-7:1e-8:1e-6,1.01e-6:1e-7:0.95e-5]'; % Diameter of dry particle (m) in ascending order
% %         %     Ddry = [1e-7:1.5e-9:1e-6,1.01e-6:1.5e-8:0.95e-5]'; % Diameter of dry particle (m) in ascending order
% %         %         Ddry = Ddry';
% %         Ddry = sortrows(Ddry,1,'ascend'); % Sort dry diameters for plotting purpose
% %         Vdry = ((3.14/6)*(Ddry.^3)).*1000000000000; % Volume of dry particle (nl)
% %
% %         Ddry = Ddry.*1000000; % Diameter of dry particle in um
% %         pd_normalDry = fitdist(log(Ddry),'Normal'); % returns mu and sigma by fitting probability distribution to data
% %
% %         %         SDdry = abs((2.302./((sqrt(2*3.14)).*(log(pd_normalDry.sigma)))).*(exp(-(((log(Ddry)-log(pd_normalDry.mu)).^2)./(2.*(log(pd_normalDry.sigma)).^2)))));
% %         SDdry = abs((2.302./((sqrt(2*3.14)).*((pd_normalDry.sigma)))).*(exp(-(((log(Ddry)-(pd_normalDry.mu)).^2)./(2.*((pd_normalDry.sigma)).^2)))));
% %
% %         DiaMeanDry = exp(pd_normalDry.mu);
% %         DiaSigmaDry = exp(pd_normalDry.sigma);

%         sheet = 2;

%  for ii = 1:11
% RH(1,ii,kkk,kkkk) = RHrange;
% lognSD = 2;n=1000;Dmu=[1.27628178
% 1.510518305
% 1.645493069
% 1.738887099
% 1.808997811
% 1.86411704
% 1.908763363
% 1.945698251
% 1.976745517
% 2.003177958
% 2.025920514
% 2.045666169
% 2.062946518
% 2.078177054
% 2.091687502
% ]';
% Ddry = DefineSizeDistribution(n,lognSD,1,Dmu(ii),0.001);

% %             prompt = 'ENTER details for size distribution --> [1,#particles,geomean of dia (um),std dev] OR [2,#particles,Min dia (um),Max dia (um)] : ';
% %             [InputSizeDistribution] = input(prompt);
% %             k=InputSizeDistribution(1);n=InputSizeDistribution(2); mu=InputSizeDistribution(3); sigma=InputSizeDistribution(4);

%                 pd_normal(:,ii,kkk,kkkk) = fitdist(log(DwetAdsIso(:,ii,kkk,kkkk)),'Normal'); % returns mu and sigma by fitting probability distribution to data
%                 DiaMean(:,ii,kkk,kkkk) = exp(pd_normal(:,ii,kkk,kkkk).mu);
%                 DiaSigma(:,ii,kkk,kkkk) = exp(pd_normal(:,ii,kkk,kkkk).sigma);
%
%                 %             SD(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*(log(pd_normal(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso(:,ii,kkk,kkkk))-log(pd_normal(:,ii,kkk,kkkk).mu)).^2)./(2.*(log(pd_normal(:,ii,kkk,kkkk).sigma)).^2)))));
%                 SD(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*((pd_normal(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso(:,ii,kkk,kkkk))-(pd_normal(:,ii,kkk,kkkk).mu)).^2)./(2.*((pd_normal(:,ii,kkk,kkkk).sigma)).^2)))));
%
%                 pd_normal1(:,ii,kkk,kkkk) = fitdist(log(DwetAdsIso1(:,ii,kkk,kkkk)),'Normal'); % returns mu and sigma by fitting probability distribution to data
%                 DiaMean1(:,ii,kkk,kkkk) = exp(pd_normal1(:,ii,kkk,kkkk).mu);
%                 DiaSigma1(:,ii,kkk,kkkk) = exp(pd_normal1(:,ii,kkk,kkkk).sigma);
%
%                 SD1(:,ii,kkk,kkkk) = abs((2.302./((sqrt(2*3.14)).*((pd_normal1(:,ii,kkk,kkkk).sigma)))).*(exp(-(((log(DwetAdsIso1(:,ii,kkk,kkkk))-(pd_normal1(:,ii,kkk,kkkk).mu)).^2)./(2.*((pd_normal1(:,ii,kkk,kkkk).sigma)).^2)))));
%
%                 pd_normalDry = fitdist(log(Ddry),'Normal'); % returns mu and sigma by fitting probability distribution to data
%                 SDdry = abs((2.302./((sqrt(2*3.14)).*((pd_normalDry.sigma)))).*(exp(-(((log(Ddry)-(pd_normalDry.mu)).^2)./(2.*((pd_normalDry.sigma)).^2)))));
%                 DiaMeanDry = exp(pd_normalDry.mu);
%                 DiaSigmaDry = exp(pd_normalDry.sigma);

% Plot the scattering angular distribution curves wrt scattering angle in polar co-ordinates
% %             figure;
% %             subplot(2,3,1);
% %             p = polarplot(angularScatPF(:,1,1,ii,kkk,kkkk),angularScatPF(:,2,1,ii,kkk,kkkk),'black');
% %             set(p,'Linewidth',1)
% %             hold on
% %             p1 = polarplot(angularScatPF1(:,1,1,ii,kkk,kkkk),angularScatPF1(:,2,1,ii,kkk,kkkk),'--r');
% %             set(p1,'Linewidth',0.75)
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(1),xSize(1,ii,kkk,kkkk),xSize1(1,ii,kkk,kkkk)));
% %             hold off
% %             subplot(2,3,2);
% %             pp = polarplot(angularScatPF(:,1,250,ii,kkk,kkkk),angularScatPF(:,2,250,ii,kkk,kkkk),'black');
% %             set(pp,'Linewidth',1)
% %             hold on
% %             pp1 = polarplot(angularScatPF1(:,1,250,ii,kkk,kkkk),angularScatPF1(:,2,250,ii,kkk,kkkk),'--r');
% %             set(pp1,'Linewidth',0.75)
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(250),xSize(250,ii,kkk,kkkk),xSize1(250,ii,kkk,kkkk)));
% %             hold off
% %             subplot(2,3,3);
% %             ppp = polarplot(angularScatPF(:,1,500,ii,kkk,kkkk),angularScatPF(:,2,500,ii,kkk,kkkk),'black');
% %             set(ppp,'Linewidth',1)
% %             hold on
% %             ppp1 = polarplot(angularScatPF1(:,1,500,ii,kkk,kkkk),angularScatPF1(:,2,500,ii,kkk,kkkk),'--r');
% %             set(ppp1,'Linewidth',0.75)
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(500),xSize(500,ii,kkk,kkkk),xSize1(500,ii,kkk,kkkk)));
% %             hold off
% %             subplot(2,3,4);
% %             ppp = polarplot(angularScatPF(:,1,750,ii,kkk,kkkk),angularScatPF(:,2,750,ii,kkk,kkkk),'black');
% %             set(ppp,'Linewidth',1)
% %             hold on
% %             ppp1 = polarplot(angularScatPF1(:,1,750,ii,kkk,kkkk),angularScatPF1(:,2,750,ii,kkk,kkkk),'--r');
% %             set(ppp1,'Linewidth',0.75)
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(750),xSize(750,ii,kkk,kkkk),xSize1(750,ii,kkk,kkkk)));
% %             hold off
% %             subplot(2,3,5);
% %             pppp = polarplot(angularScatPF(:,1,1000,ii,kkk,kkkk),angularScatPF(:,2,1000,ii,kkk,kkkk),'black');
% %             set(pppp,'Linewidth',1)
% %             hold on
% %             pppp1 = polarplot(angularScatPF1(:,1,1000,ii,kkk,kkkk),angularScatPF1(:,2,1000,ii,kkk,kkkk),'--r');
% %             set(pppp1,'Linewidth',0.75)
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(1000),xSize(1000,ii,kkk,kkkk),xSize1(1000,ii,kkk,kkkk)));
% %             hold off
% %             legend('Adsorption model','Existing models');
% %             subplot(2,3,6);
% %             txstr(1) = {sprintf('Mie Scattering Angular Distribution')};
% %             txstr(3) = {sprintf('___ Adsorption model')};
% %             txstr(4) = {sprintf('- - - Existing models (red)')};
% %             txstr(6) = {sprintf('RH = %d',RH(1,ii,kkk,kkkk))};
% %             txstr(8) = {sprintf('OIR = %d',OIR(1,1,kkk,kkkk))};
% %             txstr(10) = {sprintf('Wavelength = %d nm',lambda(1,1,1,kkkk))};
% %             text(0.1,0.6,txstr);
% %
% %             % Plot the scattering phase function wrt scattering angle
% %             figure;
% %             subplot(2,3,1);
% %             semilogy(theta,ScatPF(1,:,ii,kkk,kkkk),'black');
% %             hold on
% %             semilogy(theta,ScatPF1(1,:,ii,kkk,kkkk),'--r');
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(1),xSize(1,ii,kkk,kkkk),xSize1(1,ii,kkk,kkkk)));
% %             xlabel('Scattering angle (degrees)');
% %             ylabel('Mie scattering phase function');
% %             hold off
% %             subplot(2,3,2);
% %             semilogy(theta,ScatPF(250,:,ii,kkk,kkkk),'black');
% %             hold on
% %             semilogy(theta,ScatPF1(250,:,ii,kkk,kkkk),'--r');
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(250),xSize(250,ii,kkk,kkkk),xSize1(250,ii,kkk,kkkk)));
% %             xlabel('Scattering angle (degrees)');
% %             ylabel('Mie scattering phase function');
% %             hold off
% %             subplot(2,3,3);
% %             semilogy(theta,ScatPF(500,:,ii,kkk,kkkk),'black');
% %             hold on
% %             semilogy(theta,ScatPF1(500,:,ii,kkk,kkkk),'--r');
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(500),xSize(500,ii,kkk,kkkk),xSize1(500,ii,kkk,kkkk)));
% %             xlabel('Scattering angle (degrees)');
% %             ylabel('Mie scattering phase function');
% %             hold off
% %             subplot(2,3,4);
% %             semilogy(theta,ScatPF(750,:,ii,kkk,kkkk),'black');
% %             hold on
% %             semilogy(theta,ScatPF1(750,:,ii,kkk,kkkk),'--r');
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(750),xSize(750,ii,kkk,kkkk),xSize1(750,ii,kkk,kkkk)));
% %             xlabel('Scattering angle (degrees)');
% %             ylabel('Mie scattering phase function');
% %             hold off
% %             subplot(2,3,5);
% %             semilogy(theta,ScatPF(1000,:,ii,kkk,kkkk),'black');
% %             hold on
% %             semilogy(theta,ScatPF1(1000,:,ii,kkk,kkkk),'--r');
% %             title(sprintf('Ddry = %g um, x = [%d, %d]',Ddry(1000),xSize(1000,ii,kkk,kkkk),xSize1(1000,ii,kkk,kkkk)));
% %             xlabel('Scattering angle (degrees)');
% %             ylabel('Mie scattering phase function');
% %             hold off
% %             subplot(2,3,6);
% %             txstr(1) = {sprintf('Mie Scattering Phase Function vs Scattering angle')};
% %             txstr(3) = {sprintf('___ Adsorption model')};
% %             txstr(4) = {sprintf('- - - Existing models')};
% %             txstr(6) = {sprintf('RH = %d',RH(1,ii,kkk,kkkk))};
% %             txstr(8) = {sprintf('OIR = %d',OIR(1,1,kkk,kkkk))};
% %             txstr(10) = {sprintf('Wavelength = %d nm',lambda(1,1,1,kkkk))};
% %             text(0.05,0.6,txstr);

% %         fig1 = figure;
% %         set(fig1,'Position',[200 100 1100 650]);
% %
% %         for kkkk = 1:NN
% %             for ii = 1:N
% %
% %                 plt4 = semilogx(Ddry,SDdry,'.',DwetAdsIso1(:,ii),SD1(:,ii),'g-',DwetAdsIso(:,ii),SD(:,ii),'b--',Dwet_kappa(:,ii),SD2(:,ii),'r--'); % plot on a logarithmic x-scale
% %                 set(plt4(1),'MarkerSize',5);
% %                 set(plt4(2),'LineWidth',1.5);
% %                 set(plt4(3),'LineWidth',1.5);
% %                 set(plt4(4),'LineWidth',1.5);
% %                 hold on;
% %
% %             end
% %
% %             axis([(min(Ddry))-0.001 (max(Dwet_kappa(:,ii)))+1 0 (max(SDdry))+0.1]);
% %             title(strcat(Solute1id,'+',Solute2id,' at',RHid,'% and ',OIRid));
% %             xlabel('Particle diameter (um)');
% %             ylabel('Lognormal probability distribution function');
% %             legend('dry particles','water uptake by inorg only','water uptake by org+inorg using adsorption model','water uptake by org+inorg using kappa')
% %
% %             hold off;
% %
% %         end
% %
% %         fig2 = figure;
% %         set(fig2,'Position',[200 100 1100 650]);
% %
% %         for kkkk = 1:NN
% %             for ii = 1:N
% %
% %                 plt5 = semilogx(Ddry,gf1(:,ii),'og',Ddry,gf(:,ii),'.b',Ddry,gf_kappa(:,ii),'.r'); % plot on a logarithmic x-scale
% %                 set(plt5(1),'MarkerSize',5);
% %                 set(plt5(2),'MarkerSize',10);
% %                 set(plt5(3),'MarkerSize',10);
% %                 hold on;
% %
% %             end
% %
% %             axis([(min(Ddry))-0.001 1.4 1 4]);
% %             title(strcat(Solute1id,'+',Solute2id,' at',RHid,'% and ',OIRid));
% %             xlabel('Dry particle diameter (um)');
% %             ylabel('Growth factor');
% %             legend('water uptake by inorg only','water uptake by org+inorg using adsorption model','water uptake by org+inorg using kappa')
% %             hold off;
% %         end
% %
% %             fig3 = figure;
% %             set(fig3,'Position',[200 100 1100 650]);
% %
% %             for kkkk = 1:NN
% %                 for ii = 1:N
% %
% %                     plt6 = semilogx(aw,gf1(:,ii),'og',aw,gf(:,ii),'.b',aw,gf_kappa(:,ii),'.r'); % plot on a logarithmic x-scale
% %                     set(plt6(1),'MarkerSize',5);
% %                     set(plt6(2),'MarkerSize',10);
% %                     set(plt6(3),'MarkerSize',10);
% %                     hold on;
% %
% %                 end
% %
% %                 axis([(min(aw))-0.1 1.1 1 4]);
% %                 title(strcat(Solute1id,'+',Solute2id,' at',RHid,'% and ',OIRid));
% %                 xlabel('Water activity');
% %                 ylabel('Growth factor');
% %                 legend('water uptake by inorg only','water uptake by org+inorg using adsorption model','water uptake by org+inorg using kappa')
% %                 hold off;
% %
% %             end
% %
% %                 fig4 = figure;
% %                 set(fig4,'Position',[200 100 1100 650]);
% %
% %                 for kkkk = 1:NN
% %                     for ii = 1:N
% %
% %                         plt7 = semilogx(aw,kappa1(:,ii),'og',aw,kappa(:,ii),'.b',aw,kappa_iter(:,ii),'.r'); % plot on a logarithmic x-scale
% %                         set(plt7(1),'MarkerSize',5);
% %                         set(plt7(2),'MarkerSize',10);
% %                         set(plt7(3),'MarkerSize',10);
% %                         hold on;
% %
% %                     end
% %
% %                     axis([(min(aw))-0.1 1.1 0 1]);
% %                     title(strcat(Solute1id,'+',Solute2id,' at',RHid,'% and ',OIRid));
% %                     xlabel('Water activity');
% %                     ylabel('Hygroscopicity parameter (kappa)');
% %                     legend('water uptake by inorg only','water uptake by org+inorg using adsorption model','water uptake by org+inorg using kappa')
% %                     hold off;
% %
% %                 end
% %
% %                     % save figures
% %                     figname1 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' SizeDistribution','.fig');
% %                     saveas(fig1,figname1);
% %                     figname1 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' SizeDistribution','.jpg');
% %                     set(fig1,'PaperPositionMode','auto')
% %                     print(fig1,'-djpeg',figname1);
% %
% %                     figname2 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' DdryVSgrowthfactor','.fig');
% %                     saveas(fig2,figname2);
% %                     figname2 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' DdryVSgrowthfactor','.jpg');
% %                     set(fig2,'PaperPositionMode','auto')
% %                     print(fig2,'-djpeg',figname2);
% %
% %                     figname3 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' awVSgrowthfactor','.fig');
% %                     saveas(fig3,figname3);
% %                     figname3 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' awVSgrowthfactor','.jpg');
% %                     set(fig3,'PaperPositionMode','auto')
% %                     print(fig3,'-djpeg',figname3);
% %
% %                     figname4 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' awVSkappa','.fig');
% %                     saveas(fig4,figname4);
% %                     figname4 = strcat(folderpath,'/',Solute1id,'+', Solute2id,' awVSkappa','.jpg');
% %                     set(fig4,'PaperPositionMode','auto')
% %                     print(fig4,'-djpeg',figname4);

% %         filename1 = strcat(folderpath,'/','ScatteringCoefficient','-',Solute1id,'+', Solute2id,'-',lambdaID,'.xlsx'); % each file is for a particular wavelength, that has different OIR tabs/sheets with scat coeff versus RH
% %         writematrix(RH(:,:,kkk,kkkk),filename1,'Sheet',OIRid,'Range','A1');
% %         writematrix(TotScat(:,:,kkk,kkkk),filename1,'Sheet',OIRid,'Range','A2');
% %         writematrix(TotScat1(:,:,kkk,kkkk),filename1,'Sheet',OIRid,'Range','A3');
% %         writematrix(ScatCalcError(:,:,kkk,kkkk),filename1,'Sheet',OIRid,'Range','A4');
% %
% %         filename2 = strcat(folderpath,'/','BackScattering','-',Solute1id,'+', Solute2id,'-',lambdaID,'.xlsx'); % each file is for a particular wavelength, that has different OIR tabs/sheets with scat coeff versus RH
% %
% %         writematrix(RH(:,:,kkk,kkkk),filename2,'Sheet',OIRid,'Range','A1');
% %         writematrix(TotBackScat(:,:,kkk,kkkk),filename2,'Sheet',OIRid,'Range','A2');
% %         writematrix(TotBackScat1(:,:,kkk,kkkk),filename2,'Sheet',OIRid,'Range','A3');
% %         writematrix(BackScatCalcError(:,:,kkk,kkkk),filename2,'Sheet',OIRid,'Range','A4');
% %
% %         filename3 = strcat(folderpath,'/','AbsorptionCoefficient','-',Solute1id,'+', Solute2id,'-',lambdaID,'.xlsx'); % each file is for a particular wavelength, that has different OIR tabs/sheets with scat coeff versus RH
% %
% %         writematrix(RH(:,:,kkk,kkkk),filename3,'Sheet',OIRid,'Range','A1');
% %         writematrix(TotAbs(:,:,kkk,kkkk),filename3,'Sheet',OIRid,'Range','A2');
% %         writematrix(TotAbs1(:,:,kkk,kkkk),filename3,'Sheet',OIRid,'Range','A3');
% %         writematrix(AbsCalcError(:,:,kkk,kkkk),filename3,'Sheet',OIRid,'Range','A4');
% %
% %         filename7 = strcat(folderpath,'/','SSAlbedo','-',Solute1id,'+', Solute2id,'-',lambdaID,'.xlsx'); % each file is for a particular wavelength, that has different OIR tabs/sheets with scat coeff versus RH
% %
% %         writematrix(RH(:,:,kkk,kkkk),filename7,'Sheet',OIRid,'Range','A1');
% %         writematrix(SSA(:,:,kkk,kkkk),filename7,'Sheet',OIRid,'Range','A2');
% %         writematrix(SSA1(:,:,kkk,kkkk),filename7,'Sheet',OIRid,'Range','A3');
% %         writematrix(SSAerror(:,:,kkk,kkkk),filename7,'Sheet',OIRid,'Range','A4');
% %
% %         filename6 = strcat(folderpath,'/','MeanDiaWaterUptake','-',Solute1id,'+', Solute2id,'-','.xlsx'); % only one file for all wavelengths, that has different OIR tabs/sheets with data versus RH
% %
% %         writematrix(RH(:,:,kkk,kkkk),filename6,'Sheet',OIRid,'Range','A1');
% %         writematrix(DiaMean(:,:,kkk,kkkk),filename6,'Sheet',OIRid,'Range','A2');
% %         writematrix(DiaMean1(:,:,kkk,kkkk),filename6,'Sheet',OIRid,'Range','A3');
% %         writematrix(DiaMeanError(:,:,kkk,kkkk),filename6,'Sheet',OIRid,'Range','A4');
% %
%         filename4 = strcat(folderpath,'/','RefractiveIndex','-',Solute1id,'+', Solute2id,'-',lambdaID,'.xlsx'); % each file is for a particular wavelength, that has different OIR tabs/sheets with scat coeff versus RH
%
%         writematrix(RH(:,:,kkk,kkkk),filename4,'Sheet',OIRid,'Range','A1');
%         writematrix(RI(:,:,kkk,kkkk),filename4,'Sheet',OIRid,'Range','A2');
%         writematrix(RI1(:,:,kkk,kkkk),filename4,'Sheet',OIRid,'Range','A3');
%         writematrix(RIerror(:,:,kkk,kkkk),filename4,'Sheet',OIRid,'Range','A4');

% % TDiaSD = table(Ddry,SDdry,DwetAdsIso,SD,DwetAdsIso1,SD1);%DiaBC,
% % TDiaMean = table(RHrange,DiaMeanDry,DiaMean,DiaMean1);
% % % TRI = table(RHrange,RI,RI1,RIwaterR);
% % TsizeParam = table(xSizeDry,xSize,xSize1);
% % filename8 = strcat(folderpath,'/','DiaSD','-',Solute1id,'+', Solute2id,'-','.xlsx');
% % filename9 = strcat(folderpath,'/','DiaMean','-',Solute1id,'+', Solute2id,'-','.xlsx');
% % % filename10 = strcat(folderpath,'/','RI','-',Solute1id,'+', Solute2id,'-','.xlsx');
% % filename11 = strcat(folderpath,'/','SizeParam','-',Solute1id,'+', Solute2id,'-','.xlsx');
% % writetable(TDiaSD,filename8);
% % writetable(TDiaMean,filename9);
% % % writetable(TRI,filename10);
% % writetable(TsizeParam,filename11);

% semilogx(Ddry,SDdry);
% % DiaMean=DiaMean';
% % RI=RI';
% % TotAbs=TotAbs';