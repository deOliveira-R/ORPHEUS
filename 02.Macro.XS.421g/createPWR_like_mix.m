% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads the MICROscopic group cross sections in the Matlab
% format for materials similar to the materials of the unit cell of the
% pressurized water reactor (PWR) and calculates from them the MACROscopic
% cross sections for homogeneous mixture of dioxide uranium (fuel), natural
% mixture of zirconium isotopes (cladding) and water solution of boric acid (coolant).

function createPWR_like_mix

  
% Global structure with geometry parameters:
  global g;

  PWRmix.ng = 421;

% Path to library:
  path(path,['..' filesep '00.Lib']);
  
% Input and initialize the geometry of the PWR-like unit cell (the function
% is in '..\00.Lib')
  input_and_initialize_PWR_like;

% Path to microscopic cross section data:
  path(path,['..' filesep '01.Micro.XS.421g']);                                        % INPUT

%--------------------------------------------------------------------------  
% Call the functions for UO2 isotopes and store the data in the structures.
% As an example it is done below for temperature of 600K.
% Change when other parameters needed.
  U235 = micro_U_235__600K;                                                % INPUT
  U238 = micro_U_238__600K;                                                % INPUT
  O16 = micro_O_016__600K;                                                 % INPUT  
  PWRmix.eg = U235.eg;

% UO2 ceramic fuel is manufactured with the density lower than the
% theoretical density. The deviation is characterized with porosity which
% is the volume of voids over the total volume of the material. 0.05 (95%
% of theoretical density) is a typical value for UO2_03.
  por = 0.05;                                                              % INPUT

  
% Uranium is composed of two uranium isotopes: U235 and U238, the mass
% fraction of the U235 isotopes is called enrichment. We will used molar
% enrichment for simplicity (this is input data to be changed when needed):
  molEnrich = 0.03;                                                        % INPUT
  
% The molar fractions of U235 and U238 are as follows:
  molFrU(1) = molEnrich;
  molFrU(2) = 1 - molFrU(1);

% Mass of one "average" UO2 molecule in atomic unit mass [a.u.m.]
  UO2_03.aw = U235.aw*molFrU(1) + U238.aw*molFrU(2) + O16.aw*2.0;

% Path to material properties:
  path(path,['..' filesep '00.LIB']);
% The function returns material properties of UO2 in structure fuel
  [fuel, ~, ~] = matpro({}, {});
  
% The UO2 fuel density is theoretical density times 1 - porosity:
  UO2_03.den = fuel.rho * 1e-3 * (1 - por); % [g/cm3]
  rho = UO2_03.den*1.0e-24; % [g/(barn*cm)]
  rho = rho / 1.660538e-24; % [(a.u.m.)/(barn*cm)]
  rho = rho / UO2_03.aw; % [number of UO2 molecules/(barn*cm)]

% The names of fissionable isotopes and oxygen
  UO2_03.isoName{1} = 'U235';
  UO2_03.isoName{2} = 'U238';
  UO2_03.isoName{3} = 'O16';

% The number densities of fissionable isotopes and oxygen
  UO2_03.numDen(1) = molFrU(1) * rho * g.fuel.vFrac;
  UO2_03.numDen(2) = molFrU(2) * rho * g.fuel.vFrac;
  UO2_03.numDen(3) = 2 * rho * g.fuel.vFrac;
  
%--------------------------------------------------------------------------  
% Call the functions for Zr isotopes and store the data in the structures.
% As an example it is done below for 600K, change when other temperatures
% needed:
  Zr90 = micro_ZR090__600K;                                                % INPUT
  Zr91 = micro_ZR091__600K;                                                % INPUT
  Zr92 = micro_ZR092__600K;                                                % INPUT
  Zr94 = micro_ZR094__600K;                                                % INPUT
  Zr96 = micro_ZR096__600K;                                                % INPUT
  
  Zry.temp = 600;                                                          % INPUT
  
% Zircaloy is composed of pure Zirconium (~1% of tin neglected). There are
% four stable izotopes of zirconium: Zr090, Zr091, Zr092, Zr094 and one
% very long-lived: Zr096 with the following molar fractions:
  molFrZr = [0.5145 0.1122 0.1715 0.1738 0.0280];

% Mass of one "average" Zr atom in atomic unit mass [a.u.m.]
  Zry.aw = Zr90.aw*molFrZr(1) + Zr91.aw*molFrZr(2) + Zr92.aw*molFrZr(3) + Zr94.aw*molFrZr(4) + Zr96.aw*molFrZr(5); 

% Path to material properties:
  path(path,['..' filesep '00.LIB']);
% The function returns material properties of Zry in structure clad
  [~, ~, clad] = matpro({}, {});
  
% Zircaloy density:
  Zry.den = clad.rho*1e-3;   % [g/cm3]
  rho = Zry.den*1.0e-24;     % [g/(barn*cm)]
  rho = rho / 1.660538e-24;  % [(a.u.m.)/(barn*cm)]
  rho = rho / Zry.aw;        % [number of Zr atoms/(barn*cm)]

% The names of isotopes
  Zry.isoName{1} = 'Zr90';
  Zry.isoName{2} = 'Zr91';
  Zry.isoName{3} = 'Zr92';
  Zry.isoName{4} = 'Zr94';
  Zry.isoName{5} = 'Zr96';

% The number densities of isotopes:
  Zry.numDen(1) = molFrZr(1) * rho * g.clad.vFrac;
  Zry.numDen(2) = molFrZr(2) * rho * g.clad.vFrac;
  Zry.numDen(3) = molFrZr(3) * rho * g.clad.vFrac;
  Zry.numDen(4) = molFrZr(4) * rho * g.clad.vFrac;
  Zry.numDen(5) = molFrZr(5) * rho * g.clad.vFrac;
  
%--------------------------------------------------------------------------  
% Call the functions for H2O and B isotopes and store the data in the
% structures. As an example it is done below for temperature of 600K, 
% pressure of 16 MPa and boron concentration of 4000 ppm.
% Change when other parameters needed.
  H01 = micro_H_001__600K;                                                 % INPUT
  O16_ = micro_O_016__600K;                                                % INPUT
  B10 = micro_B_010__600K;                                                 % INPUT
  B11 = micro_B_011__600K;                                                 % INPUT
  
  H2OB.temp = 600; %K                                                      % INPUT
  H2OB.p = 16; %MPa                                                        % INPUT
  H2OB.bConc = 4000e-6; % 1e-6 = 1 ppm                                     % INPUT
  
% Boron is composed of two stable isotopes: B10 and B11 with the following
% molar fractions:
  molFrB = [0.199 0.801];

% Mass of one "average" H2OB molecule in atomic unit mass [a.u.m.]:
  H2OB.aw = 2*H01.aw + O16_.aw + H2OB.bConc * (molFrB(1)*B10.aw + molFrB(2)*B11.aw);

% Path to steam-water properties:
  path(path,['..' filesep '00.XSteam']);                                               % INPUT
% The function returns water density at specified pressure (MPa) and 
% temperature (C):
  density = XSteam('rho_pt', H2OB.p*10, H2OB.temp-273);

% The water density:
  H2OB.den = density *1e-3;  % [g/cm3]
  rho = H2OB.den*1.0e-24;    % [g/(barn*cm)]
  rho = rho / 1.660538e-24;  % [(a.u.m.)/(barn*cm)]
  rho = rho / H2OB.aw;       % [number of H2O molecules/(barn*cm)]

% The names of isotopes
  H2OB.isoName{1} = 'H01';
  H2OB.isoName{2} = 'O16';
  H2OB.isoName{3} = 'B10';
  H2OB.isoName{4} = 'B11';

% The number densities of isotopes:
  H2OB.numDen(1) = 2 * rho * g.cool.vFrac;
  H2OB.numDen(2) = rho * g.cool.vFrac;
  H2OB.numDen(3) = rho*H2OB.bConc*molFrB(1) * g.cool.vFrac;
  H2OB.numDen(4) = rho*H2OB.bConc*molFrB(2) * g.cool.vFrac;

%--------------------------------------------------------------------------  
% Prepare for sigma-zero iterations:
  sigTtab = {U235.sigT, U238.sigT, O16.sigT, Zr90.sigT, Zr91.sigT, Zr92.sigT, Zr94.sigT, Zr96.sigT, H01.sigT, O16_.sigT, B10.sigT, B11.sigT};
  sig0tab = {U235.sig0, U238.sig0, O16.sig0, Zr90.sig0, Zr91.sig0, Zr92.sig0, Zr94.sig0, Zr96.sig0, H01.sig0, O16_.sig0, B10.sig0, B11.sig0};
  aDen = [UO2_03.numDen'; Zry.numDen'; H2OB.numDen'];
  
% SigEscape -- escape cross section, for simple convex objects (such as
% plates, spheres, or cylinders) is given by S/(4V), where V and S are the
% volume and surface area of the object, respectively
  SigEscape = 0;

  fprintf('Sigma-zero iterations. ');
  [PWRmix.sig0] = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape);
  fprintf('Done.\n');
  
  fprintf('Interpolation of microscopic cross sections for the found sigma-zeros. ');
  sigCtab = {U235.sigC, U238.sigC, O16.sigC, Zr90.sigC, Zr91.sigC, Zr92.sigC, Zr94.sigC, Zr96.sigC, H01.sigC, O16_.sigC, B10.sigC, B11.sigC};
  sigLtab = {U235.sigL, U238.sigL, O16.sigL, Zr90.sigL, Zr91.sigL, Zr92.sigL, Zr94.sigL, Zr96.sigL, H01.sigL, O16_.sigL, B10.sigL, B11.sigL};
  sigFtab = {U235.sigF, U238.sigF, O16.sigF, Zr90.sigF, Zr91.sigF, Zr92.sigF, Zr94.sigF, Zr96.sigF, H01.sigF, O16_.sigF, B10.sigF, B11.sigF};
  sig2 =    {U235.sig2, U238.sig2, O16.sig2, Zr90.sig2, Zr91.sig2, Zr92.sig2, Zr94.sig2, Zr96.sig2, H01.sig2, O16_.sig2, B10.sig2, B11.sig2};
  for ig = 1:PWRmix.ng
    % Number of isotopes in the mixture
      nIso = length(aDen);
    % Loop over isotopes
      for iIso = 1:nIso
        % Find cross sections for the found sigma-zeros
          if length(sig0tab{iIso}) == 1
             sigC(iIso,ig) = sigCtab{iIso}(1,ig);
             sigL(iIso,ig) = sigLtab{iIso}(1,ig);
             sigF(iIso,ig) = sigFtab{iIso}(1,ig);
          else
             log10sig0 = min(10,max(0,log10(PWRmix.sig0(iIso,ig))));
             sigC(iIso,ig) = interp1(log10(sig0tab{iIso})', sigCtab{iIso}(:,ig), log10sig0);
             sigL(iIso,ig) = interp1(log10(sig0tab{iIso})', sigLtab{iIso}(:,ig), log10sig0);
             sigF(iIso,ig) = interp1(log10(sig0tab{iIso})', sigFtab{iIso}(:,ig), log10sig0);
          end
      end
  end
 
 % Find scattering matrices for the found sigma-zeros
  for j = 0:2
      sigS{1+j,1}  = interpSigS(1+j, U235, PWRmix.sig0(1,:));
      sigS{1+j,2}  = interpSigS(1+j, U238, PWRmix.sig0(2,:));
      sigS{1+j,3}  = interpSigS(1+j, O16,  PWRmix.sig0(3,:));
      sigS{1+j,4}  = interpSigS(1+j, Zr90, PWRmix.sig0(4,:));
      sigS{1+j,5}  = interpSigS(1+j, Zr91, PWRmix.sig0(5,:));
      sigS{1+j,6}  = interpSigS(1+j, Zr92, PWRmix.sig0(6,:));
      sigS{1+j,7}  = interpSigS(1+j, Zr94, PWRmix.sig0(7,:));
      sigS{1+j,8}  = interpSigS(1+j, Zr96, PWRmix.sig0(8,:));
      sigS{1+j,9}  = interpSigS(1+j, H01,  PWRmix.sig0(9,:));
      sigS{1+j,10} = interpSigS(1+j, O16_, PWRmix.sig0(10,:));
      sigS{1+j,11} = interpSigS(1+j, B10,  PWRmix.sig0(11,:));
      sigS{1+j,12} = interpSigS(1+j, B11,  PWRmix.sig0(12,:));
  end
  fprintf('Done.\n');
 
%--------------------------------------------------------------------------  
% Macroscopic cross section [1/cm] is microscopic cross section for the 
% "average" molecule [barn] times the number density [number of
% molecules/(barn*cm)]
  PWRmix.SigC = sigC'*aDen;
  PWRmix.SigL = sigL'*aDen;
  PWRmix.SigF = sigF'*aDen;
  PWRmix.SigP = U235.nubar.*sigF(1,:)*aDen(1) + U238.nubar.*sigF(2,:)*aDen(2);
  for j = 0:2
      PWRmix.SigS{1+j} = sigS{1+j,1}*aDen(1) + sigS{1+j,2}*aDen(2) + sigS{1+j,3}*aDen(3) + sigS{1+j,4}*aDen(4) + sigS{1+j,5}*aDen(5) + sigS{1+j,6}*aDen(6) + sigS{1+j,7}*aDen(7) + sigS{1+j,8}*aDen(8) + sigS{1+j,9}*aDen(9) + sigS{1+j,10}*aDen(10) + sigS{1+j,11}*aDen(11) + sigS{1+j,12}*aDen(12);
  end
  PWRmix.Sig2 = sig2{1}*aDen(1) + sig2{2}*aDen(2) + sig2{3}*aDen(3) + sig2{4}*aDen(4) + sig2{5}*aDen(5) + sig2{6}*aDen(6) + sig2{7}*aDen(7) + sig2{8}*aDen(8) + sig2{9}*aDen(9) + sig2{10}*aDen(10) + sig2{11}*aDen(11) + sig2{12}*aDen(12);
  PWRmix.SigT = PWRmix.SigC + PWRmix.SigL + PWRmix.SigF + sum(PWRmix.SigS{1+0},2) + sum(PWRmix.Sig2,2);

% For simplicity fission spectrum of the mixture assumed equal to fission spectrum of U235
  PWRmix.chi =  U235.chi;
% Number of fissile isotopes, macroscopic production cross section and fission spectrum for every fissile isotopes
  PWRmix.nFis = 2;
  PWRmix.fis.SigP(1,:) = U235.nubar.*sigF(1,:)*aDen(1);
  PWRmix.fis.SigP(2,:) = U238.nubar.*sigF(2,:)*aDen(2);
  PWRmix.fis.chi(1,:) = U235.chi;
  PWRmix.fis.chi(2,:) = U238.chi;
  
% Make a file name which include the isotope name
  matName = 'macro421_PWR_like_mix';
  
% Make a header for the file to be created with important parameters for
% which the macroscopic cross sections were generated
  PWRmix.header{1} = sprintf('%% ---------------------------------------------------------');
  PWRmix.header{2} = sprintf('%% Matlab-based Open-source Reactor Physics Education System');
  PWRmix.header{3} = sprintf('%% ---------------------------------------------------------');
  PWRmix.header{4} = sprintf('%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.');
  PWRmix.header{5} = sprintf('%%');
  PWRmix.header{6} = sprintf('%% Macroscopic cross sections for homogeneos mixture of PWR unit cell materials');
 
% Just write zeros for aw, den and temp
  PWRmix.aw = 0;
  PWRmix.den = 0;
  PWRmix.temp = 0;
  
  PWRmix.numDen = aDen;
  PWRmix.isoName = [UO2_03.isoName'; Zry.isoName'; H2OB.isoName'];

% Change the units of number density from 1/(barn*cm) to 1/cm2
  PWRmix.numDen = PWRmix.numDen*1e24;

% Finally create the file with macroscopic cross sections
  writeMacroXS(PWRmix,matName);
  
  fprintf("\nU235 %f \n", UO2_03.numDen(1))
  fprintf("U238 %f \n", UO2_03.numDen(2))
  fprintf("O16 %f \n", UO2_03.numDen(3))
  
  fprintf("\nZr90 %f \n", Zry.numDen(1))
  fprintf("Zr91 %f \n", Zry.numDen(2))
  fprintf("Zr92 %f \n", Zry.numDen(3))
  fprintf("Zr94 %f \n", Zry.numDen(4))
  fprintf("Zr96 %f \n", Zry.numDen(5))
  
  fprintf("\nH2 %f\n", H2OB.numDen(1))
  fprintf("O16 %f\n", H2OB.numDen(2))
  fprintf("B10 %f\n", H2OB.numDen(3))
  fprintf("B11 %f\n", H2OB.numDen(4))
end