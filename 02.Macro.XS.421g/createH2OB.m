% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads the MICROscopic group cross sections in the Matlab
% format and calculates from them the MACROscopic cross sections for water
% solution of boric acid which is similar to the coolant of the pressurized
% water reactor.

function createH2OB

% number of energy groups
  H2OB.ng = 421;
  
% Boron is composed of two stable isotopes: B10 and B11 with the following
% molar fractions:
  molFrB = [0.199 0.801];
  
% Path to library:
  path(path,['..' filesep '00.Lib']);
  
% Input and initialize the geometry of the PWR-like unit cell (the function
% is in '..\00.Lib')
  input_and_initialize_PWR_like;

% Path to microscopic cross section data:
  path(path,['..' filesep '01.Micro.XS.421g']);

% Call the functions for H2O and B isotopes and store the data in the
% structures. As an example it is done below for temperature of 600K, 
% pressure of 16 MPa and boron concentration of 4000 ppm.
% Change when other parameters needed.
  H01 = micro_H_001__600K;                                                 % INPUT
  O16 = micro_O_016__600K;                                                 % INPUT
  B10 = micro_B_010__600K;                                                 % INPUT
  B11 = micro_B_011__600K;                                                 % INPUT
  
  H2OB.temp = 600; %K                                                      % INPUT
  H2OB.p = 16; %MPa                                                        % INPUT
  H2OB.bConc = 4000e-6; % 1e-6 = 1 ppm                                     % INPUT
  H2OB.eg = H01.eg; 

% Mass of one "average" H2OB molecule in atomic unit mass [a.u.m.]:
  H2OB.aw = 2*H01.aw + O16.aw + H2OB.bConc * (molFrB(1)*B10.aw + molFrB(2)*B11.aw);

% Path to steam-water properties:
  path(path,['..' filesep '00.XSteam']);
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
  H2OB.numDen(1) = 2 * rho;
  H2OB.numDen(2) = rho;
  H2OB.numDen(3) = rho*H2OB.bConc*molFrB(1);
  H2OB.numDen(4) = rho*H2OB.bConc*molFrB(2);

% Prepare for sigma-zero iterations:
  sigTtab = {H01.sigT, O16.sigT, B10.sigT, B11.sigT};
  sig0tab = {H01.sig0, O16.sig0, B10.sig0, B11.sig0};
  aDen = H2OB.numDen';

% SigEscape -- escape cross section, for simple convex objects (such as
% plates, spheres, or cylinders) is given by S/(4V), where V and S are the
% volume and surface area of the object, respectively
  SigEscape = 0;

  fprintf('Sigma-zero iterations. ');
  H2OB.sig0 = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape);
  fprintf('Done.\n');

  fprintf('Interpolation of microscopic cross sections for the found sigma-zeros. ');
  sigCtab = {H01.sigC, O16.sigC, B10.sigC, B11.sigC};
  sigLtab = {H01.sigL, O16.sigL, B10.sigL, B11.sigL};
  for ig = 1:H2OB.ng
    % Loop over isotopes
      for iIso = 1:4
        % Find cross sections for the sigma-zero
          if length(sig0tab{iIso}) == 1
             sigC(iIso,ig) = sigCtab{iIso}(1,ig);
             sigL(iIso,ig) = sigLtab{iIso}(1,ig);
          else
             log10sig0 = min(10,max(0,log10(H2OB.sig0(iIso,ig))));
             sigC(iIso,ig) = interp1(log10(sig0tab{iIso})', sigCtab{iIso}(:,ig), log10sig0);
             sigL(iIso,ig) = interp1(log10(sig0tab{iIso})', sigLtab{iIso}(:,ig), log10sig0);
          end
      end
  end
  
  for j = 0:2
      sigS{1+j,1} = interpSigS(1+j, H01, H2OB.sig0(1,:));
      sigS{1+j,2} = interpSigS(1+j, O16, H2OB.sig0(2,:));
      sigS{1+j,3} = interpSigS(1+j, B10, H2OB.sig0(3,:));
      sigS{1+j,4} = interpSigS(1+j, B11, H2OB.sig0(4,:));
  end
  fprintf('Done.\n');
 
% Macroscopic cross section [1/cm] is microscopic cross section for the 
% molecule [barn] times the number density [number of molecules/(barn*cm)]
  H2OB.SigC = sigC'*aDen;
  H2OB.SigL = sigL'*aDen; 
  for j = 0:2
      H2OB.SigS{1+j} = sigS{1+j,1}*aDen(1) + sigS{1+j,2}*aDen(2) + sigS{1+j,3}*aDen(3) + sigS{1+j,4}*aDen(4); 
  end
  H2OB.Sig2 = H01.sig2*aDen(1) + O16.sig2*aDen(2) + B10.sig2*aDen(3) + B11.sig2*aDen(4);
  H2OB.SigT = H2OB.SigC + H2OB.SigL + sum(H2OB.SigS{1+0},2) + sum(H2OB.Sig2,2);
  
  H2OB.SigP(1) = 0.0;

% Make a file name which include the isotope name and the temperature
  if H2OB.temp < 1000
     matName = ['macro421_H2OB__',num2str(round(H2OB.temp)),'K']; % name of the file with a temperature index
  else
     matName = ['macro421_H2OB_',num2str(round(H2OB.temp)),'K']; % name of the file with a temperature index
  end
  
% Make a header for the file to be created with important parameters for
% which the macroscopic cross sections were generated
  H2OB.header{1} = sprintf('%% ---------------------------------------------------------');
  H2OB.header{2} = sprintf('%% Matlab-based Open-source Reactor Physics Education System');
  H2OB.header{3} = sprintf('%% ---------------------------------------------------------');
  H2OB.header{4} = sprintf('%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.');
  H2OB.header{5} = sprintf('%%');
  H2OB.header{6} = sprintf('%% Macroscopic cross sections for water solution of boric acid');
  H2OB.header{7} = sprintf('%% Water temperature:   %6.1f K',H2OB.temp);
  H2OB.header{8} = sprintf('%% Water pressure:      %4.1f MPa',H2OB.p);
  H2OB.header{9} = sprintf('%% Water density:       %8.5f g/cm3',H2OB.den);
  H2OB.header{10} = sprintf('%% Boron concentration: %6.1f ppm',H2OB.bConc*1e6);
  
% Change the units of number density from 1/(barn*cm) to 1/cm2
  H2OB.numDen = H2OB.numDen*1e24;

% Finally create the file with macroscopic cross sections
  writeMacroXS(H2OB,matName);
 end