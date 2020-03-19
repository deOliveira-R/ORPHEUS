% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads the MICROscopic group cross sections in the Matlab
% format and calculates from them the MACROscopic cross sections for water
% solution of uranium-235 which is a homogeneous aqueous reactor.

function createH2OU

% number of energy groups
  H2OU.ng = 421;
  
% Path to microscopic cross section data:
  path(path,['..' filesep '01.Micro.XS.421g']);                                        % INPUT

% Call the functions for H2O and B isotopes and store the data in the
% structures. As an example it is done below for temperature of 294K, 
% pressure of 7 MPa and boron concentration of 760 ppm.
% Change when other parameters needed.
  H01 = micro_H_001__294K;                                                 % INPUT
  O16 = micro_O_016__294K;                                                 % INPUT
  U235 = micro_U_235__294K;                                                % INPUT
  
  H2OU.temp = 294; %K                                                      % INPUT
  H2OU.p = 0.1; %MPa                                                       % INPUT
  H2OU.Uconc = 1000e-6; % 1e-6 = 1 ppm                                     % INPUT
  H2OU.eg = H01.eg;

% Mass of one "average" H2OU molecule in atomic unit mass [a.u.m.]:
  H2OU.aw = 2*H01.aw + O16.aw + H2OU.Uconc*U235.aw;

% Path to steam-water properties:
  path(path,['..' filesep '00.XSteam']);                                               % INPUT
% The function returns water density at specified pressure (MPa) and 
% temperature (C):
  density = XSteam('rho_pt', H2OU.p/10, H2OU.temp-273);

% The water density:
  H2OU.den = density *1e-3;  % [g/cm3]
  rho = H2OU.den*1.0e-24;    % [g/(barn*cm)]
  rho = rho / 1.660538e-24;  % [(a.u.m.)/(barn*cm)]
  rho = rho / H2OU.aw;       % [number of H2O molecules/(barn*cm)]

% The names of fissionable isotopes and oxygen
  H2OU.isoName{1} = 'H01';
  H2OU.isoName{2} = 'O16';
  H2OU.isoName{3} = 'U235';

% The number densities of isotopes:
  H2OU.numDen(1) = 2 * rho;
  H2OU.numDen(2) = rho;
  H2OU.numDen(3) = rho*H2OU.Uconc;

% Prepare for sigma-zero iterations:
  sigTtab = {H01.sigT, O16.sigT, U235.sigT};
  sig0tab = {H01.sig0, O16.sig0, U235.sig0};
  aDen = H2OU.numDen';

% SigEscape -- escape cross section, for simple convex objects (such as
% plates, spheres, or cylinders) is given by S/(4V), where V and S are the
% volume and surface area of the object, respectively
  SigEscape = 0;

  fprintf('Sigma-zero iterations. ');
  H2OU.sig0 = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape);
  fprintf('Done.\n');

  fprintf('Interpolation of microscopic cross sections for the found sigma-zeros. ');
  sigCtab = {H01.sigC, O16.sigC, U235.sigC};
  sigLtab = {H01.sigL, O16.sigL, U235.sigL};
  sigFtab = {H01.sigF, O16.sigF, U235.sigF};
  for ig = 1:H2OU.ng
    % Loop over isotopes
      for iIso = 1:3
        % Find cross sections for the sigma-zero
          if length(sig0tab{iIso}) == 1
             sigC(iIso,ig) = sigCtab{iIso}(1,ig);
             sigL(iIso,ig) = sigLtab{iIso}(1,ig);
             sigF(iIso,ig) = sigFtab{iIso}(1,ig);
          else
             log10sig0 = min(10,max(0,log10(H2OU.sig0(iIso,ig))));
             sigC(iIso,ig) = interp1(log10(sig0tab{iIso})', sigCtab{iIso}(:,ig), log10sig0);
             sigL(iIso,ig) = interp1(log10(sig0tab{iIso})', sigLtab{iIso}(:,ig), log10sig0);
             sigF(iIso,ig) = interp1(log10(sig0tab{iIso})', sigFtab{iIso}(:,ig), log10sig0);
          end
      end
  end
  
  for j = 0:2
      sigS{1+j,1} = interpSigS(1+j, H01, H2OU.sig0(1,:));
      sigS{1+j,2} = interpSigS(1+j, O16, H2OU.sig0(2,:));
      sigS{1+j,3} = interpSigS(1+j, U235, H2OU.sig0(3,:));
  end
  fprintf('Done.\n');
 
% Macroscopic cross section [1/cm] is microscopic cross section for the 
% molecule [barn] times the number density [number of molecules/(barn*cm)]
  H2OU.SigC = sigC'*aDen;
  H2OU.SigL = sigL'*aDen;
  H2OU.SigF = sigF'*aDen;
  H2OU.SigP = U235.nubar.*sigF(3,:)*aDen(3);
  for j = 0:2
      H2OU.SigS{1+j} = sigS{1+j,1}*aDen(1) + sigS{1+j,2}*aDen(2) + sigS{1+j,3}*aDen(3); 
  end
  H2OU.Sig2 = H01.sig2*aDen(1) + O16.sig2*aDen(2) + U235.sig2*aDen(3); 
  H2OU.SigT = H2OU.SigC + H2OU.SigL + H2OU.SigF + sum(H2OU.SigS{1+0},2) + sum(H2OU.Sig2,2);
  
% Fission spectrum
  H2OU.chi = U235.chi;
  
% Make a file name which include the isotope name and the temperature
  if H2OU.temp < 1000
     matName = ['macro421_H2OU__',num2str(round(H2OU.temp)),'K']; % name of the file with a temperature index
  else
     matName = ['macro421_H2OU_',num2str(round(H2OU.temp)),'K']; % name of the file with a temperature index
  end
  
% Make a header for the file to be created with important parameters for
% which the macroscopic cross sections were generated
  H2OU.header{1} = sprintf('%% ---------------------------------------------------------');
  H2OU.header{2} = sprintf('%% Matlab-based Open-source Reactor Physics Education System');
  H2OU.header{3} = sprintf('%% ---------------------------------------------------------');
  H2OU.header{4} = sprintf('%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.');
  H2OU.header{5} = sprintf('%%');
  H2OU.header{6} = sprintf('%% Macroscopic cross sections for water solution of uranium-235');
  H2OU.header{7} = sprintf('%% Water temperature:   %6.1f K',H2OU.temp);
  H2OU.header{8} = sprintf('%% Water pressure:      %4.1f MPa',H2OU.p);
  H2OU.header{9} = sprintf('%% Water density:       %8.5f g/cm3',H2OU.den);
  H2OU.header{10} = sprintf('%% U-235 concentration:  %6.1f ppm',H2OU.Uconc*1e6);
  
% Change the units of number density from 1/(barn*cm) to 1/cm2
  H2OU.numDen = H2OU.numDen*1e24;

% Finally create the file with macroscopic cross sections
  writeMacroXS(H2OU,matName);
 end