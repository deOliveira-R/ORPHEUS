% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads the MICROscopic group cross sections in the Matlab
% format and calculates from them the MACROscopic cross sections for
% dioxide uranium which is the fuel material of the pressurized water 
% reactor.

function createUO2_03

% Global structure with geometry parameters
  global g

% number of energy groups
  UO2_03.ng = 421;
  
% Path to library:
  path(path,['..' filesep '00.Lib']);
  
% input and initialize the geometry of the PWR unit cell (the function is
% in '..\00.Lib')
  input_and_initialize_PWR_like;

% Path to microscopic cross section data:
  path(path,['..' filesep '01.Micro.XS.421g']);                                        % INPUT
  
% Call the functions for UO2 isotopes and store the data in the structures.
% As an example it is done below for temperature of 900 K.
% Change when other parameters needed.
  U235 = micro_U_235__900K;                                                % INPUT
  U238 = micro_U_238__900K;                                                % INPUT
  O16 = micro_O_016__900K;                                                 % INPUT
  
  UO2_03.temp = 900;                                                       % INPUT
  UO2_03.eg = U235.eg;

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
  rho = UO2_03.den*1.0e-24;                 % [g/(barn*cm)]
  rho = rho / 1.660538e-24;                 % [(a.u.m.)/(barn*cm)]
  rho = rho / UO2_03.aw;                    % [number of UO2 molecules/(barn*cm)]

% The names of fissionable isotopes and oxygen
  UO2_03.isoName{1} = 'U235';
  UO2_03.isoName{2} = 'U238';
  UO2_03.isoName{3} = 'O16';

% The number densities of fissionable isotopes and oxygen
  UO2_03.numDen(1) = molFrU(1) * rho;
  UO2_03.numDen(2) = molFrU(2) * rho;
  UO2_03.numDen(3) = 2 * rho;
  
% Prepare for sigma-zero iterations:
  sigTtab = {U235.sigT, U238.sigT, O16.sigT};
  sig0tab = {U235.sig0, U238.sig0, O16.sig0};
  aDen = UO2_03.numDen';
  
% SigEscape -- escape cross section, for simple convex objects (such as
% plates, spheres, or cylinders) is given by S/(4V), where V and S are the
% volume and surface area of the object, respectively
  SigEscape = 1 / (2*g.fuel.rOut*100);

  fprintf('Sigma-zero iterations. ');
  UO2_03.sig0 = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape);
  fprintf('Done.\n');
  
  fprintf('Interpolation of microscopic cross sections for the found sigma-zeros. ');
  sigCtab = {U235.sigC, U238.sigC, O16.sigC};
  sigLtab = {U235.sigL, U238.sigL, O16.sigL};
  sigFtab = {U235.sigF, U238.sigF, O16.sigF};
  for ig = 1:UO2_03.ng
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
             log10sig0 = min(10,max(0,log10(UO2_03.sig0(iIso,ig))));
             sigC(iIso,ig) = interp1(log10(sig0tab{iIso})', sigCtab{iIso}(:,ig), log10sig0);
             sigL(iIso,ig) = interp1(log10(sig0tab{iIso})', sigLtab{iIso}(:,ig), log10sig0);
             sigF(iIso,ig) = interp1(log10(sig0tab{iIso})', sigFtab{iIso}(:,ig), log10sig0);
          end
      end
  end
 
 % Find scattering matrices for the sigma-zeros
  for j = 0:1
      sigS{1+j,1} = interpSigS(1+j, U235, UO2_03.sig0(1,:));
      sigS{1+j,2} = interpSigS(1+j, U238, UO2_03.sig0(2,:));
      sigS{1+j,3} = interpSigS(1+j, O16, UO2_03.sig0(3,:));
  end
  fprintf('Done.\n');
 
% Macroscopic cross section [1/cm] is microscopic cross section for the 
% "average" molecule [barn] times the number density [number of
% molecules/(barn*cm)]
  UO2_03.SigC = sigC'*aDen;
  UO2_03.SigL = sigL'*aDen;
  UO2_03.SigF = sigF'*aDen;
  UO2_03.SigP = U235.nubar.*sigF(1,:)*aDen(1) + U238.nubar.*sigF(2,:)*aDen(2);
  for j = 0:1
      UO2_03.SigS{1+j} = sigS{1+j,1}*aDen(1) + sigS{1+j,2}*aDen(2) + sigS{1+j,3}*aDen(3);
  end
  UO2_03.Sig2 = U235.sig2*aDen(1) + U238.sig2*aDen(2) + O16.sig2*aDen(3);
  UO2_03.SigT = UO2_03.SigC + UO2_03.SigL + UO2_03.SigF + sum(UO2_03.SigS{1+0},2) + sum(UO2_03.Sig2,2);

% For simplicity assume fission spectrum of the mixture equal to fission
% spectrum of U235
  UO2_03.chi =  U235.chi;
    
% Make a file name which include the isotope name and the temperature
  if UO2_03.temp < 1000
     matName = ['macro421_UO2_03__',num2str(round(UO2_03.temp)),'K']; % name of the file with a temperature index
  else
     matName = ['macro421_UO2_03_',num2str(round(UO2_03.temp)),'K']; % name of the file with a temperature index
  end
  
% Make a header for the file to be created with important parameters for
% which the macroscopic cross sections were generated
  UO2_03.header{1} = sprintf('%% ---------------------------------------------------------');
  UO2_03.header{2} = sprintf('%% Matlab-based Open-source Reactor Physics Education System');
  UO2_03.header{3} = sprintf('%% ---------------------------------------------------------');
  UO2_03.header{4} = sprintf('%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.');
  UO2_03.header{5} = sprintf('%%');
  UO2_03.header{6} = sprintf('%% Macroscopic cross sections for uranium dioxide');
  UO2_03.header{7} = sprintf('%% UO2 temperature:    %6.1f K',UO2_03.temp);
  UO2_03.header{8} = sprintf('%% UO2 porosity:       %4.1f %%',por*100);
  UO2_03.header{9} = sprintf('%% UO2 density:        %7.3f g/cm3',UO2_03.den);
  UO2_03.header{10} = sprintf('%% Enrichment by U235: %4.1f mol%%',molEnrich*100);

% Change the units of number density from 1/(barn*cm) to 1/cm2
  UO2_03.numDen = UO2_03.numDen*1e24;

% Finally create the file with macroscopic cross sections
  writeMacroXS(UO2_03,matName);
end