% --------------------------------
% Reactor Physics MATLAB Functions
% --------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function specifies some material properties of uranium dioxide,
% inert gases and zircaloy based on openly available documents.

function [fuel, gap, clad] = matpro(fuel, gap, clad)

% UO2: theoretical density (kg/m3) MATPRO(2003) p. 2-56
  fuel.rho = 10980; 
% UO2: specific heat (J/kg-K)
  fuel.cp = @(T) 162.3 + 0.3038.*T - 2.391e-4.*T.^2 + 6.404e-8.*T.^3;
% UO2 thermal conductivity (W/m-K) MATPRO
  fuel.k = @(T, Bu, por) ( 1 ./ (0.0452 + 0.000246*T + 0.00187*Bu   + 0.038*(1 - 0.9*exp(-0.04*Bu)).*Bu.^0.28./(1 + 396*exp(-6380./T))) + 3.5e9*exp(-16360./T)./T.^2 ) * 1.0789 .* (1 - por)./(1 + por/2);
% UO2: thermal expansion (m/m) MATPRO(2003)
  fuel.thExp = @(T) (T/1000 - 0.3 + 4*exp(-5000 ./ T))/100;
% UO2: Young's modulus (MPa) MATPRO(2003) p. 2-58
  fuel.E = @(T,por) 2.334e5 * (1 - 2.752*por).*(1 - 1.0915e-4*T);
% UO2: Poisson ratio (-) MATPRO(2003) p. 2-68
  fuel.nu = 0.316;
% UO2: swelling rate (1/s) MATPRO(2003)
  fuel.swelRate = @(dFdt,F,T) 2.5e-29*dFdt + (T < 2800).*(8.8e-56 * dFdt .* (2800 - T).^11.73 .* exp(-0.0162*(2800 - T)) .* exp(-8e-27 * F));
% UO2: thermal creep rate (1/s) a simplified correlation for sensitivity
% study
  fuel.creepRate = @(sig,T) 5e5 *sig .* exp( -4e5./(8.314*T) );

%--------------------------------------------------------------------------
% He: thermal conductivity (W/m-K)
  gap.kHe = @(T) 2.639e-3*T.^0.7085;
% Xe: thermal conductivity (W/m-K)
  gap.kXe = @(T) 4.351e-5*T.^0.8616;
% Kr: thermal conductivity (W/m-K)
  gap.kKr = @(T) 8.247e-5*T.^0.8363; 
% auxiliary function for gas mixture gas conductivity calculation (-) MATPRO
%  psi = @(k1,k2,M1,M2) (1 + sqrt( sqrt(M1/M2)*k1./k2 )).^2  * ( 1 + 2.41*(M1-M2)*(M1-0.142*M2)/sqrt(M1+M2) ) / sqrt(8*(1+M1/M2));
  psi = @(k1,k2,M1,M2) (1 + sqrt( sqrt(M1/M2)*k1./k2 )).^2  / sqrt(8*(1+M1/M2));
% gas mixture gas conductivity (-) MATPRO
  gap.kGasMixFun = @(k,x,M) k{1}*x(1) ./ ( psi(k{1},k{1},M(1),M(1))*x(1) + psi(k{1},k{2},M(1),M(2))*x(2) + psi(k{1},k{3},M(1),M(3))*x(3) ) + ... % gas mixture gas conductivity (W/m-K) MATPRO
                            k{2}*x(2) ./ ( psi(k{2},k{1},M(2),M(1))*x(1) + psi(k{2},k{2},M(2),M(2))*x(2) + psi(k{2},k{3},M(2),M(3))*x(3) ) + ...
                            k{3}*x(3) ./ ( psi(k{3},k{1},M(3),M(1))*x(1) + psi(k{3},k{2},M(3),M(2))*x(2) + psi(k{3},k{3},M(3),M(3))*x(3) );
%--------------------------------------------------------------------------
% Zry: density (kg/m3)
  clad.rho = 6600;
% Zry: specific heat (J/kg-K)
  clad.cp = @(T) 252.54 + 0.11474.*T;
% Zry thermal conductivity (W/m-K)
  clad.k = @(T) 7.51 + 2.09e-2.*T - 1.45e-5.*T.^2 + 7.67e-9.*T.^3;
% Zry: thermal expansion (-) PNNL(2010) p.3-16
  clad.thExp = @(T) -2.373e-5 + (T-273.15)*6.721e-6;
% Zry: Young's modulus (MPa) PNNL(2010) p. 3-20 (cold work assumed zero)
  clad.E = @(T) 1.088e5 - 54.75*T;
% Zry: Poisson ratio (-) MATPRO(2003) p. 4-242
  clad.nu = 0.3;
% Zry: thermal creep rate (1/s) a simplified correlation for sensitivity
% study
  clad.creepRate = @(sig,T) 1e5 *sig .* exp( -2e5./(8.314*T) );
% Zry: strength coefficient MATPRO
  clad.K = @(T) (T<743).*(2.257e9 + T.*(-5.644e6 + T.*(7.525e3 - T*4.33167))) + (T>=743 & T<1090).*(2.522488e6*exp(2.8500027e6./T.^2)) + (T>=1090 & T<1255).*(184.1376039e6 - 1.4345448e5*T) + (T>=1255).*(4.330e7 + T.*(-6.685e4 + T.*(37.579 - T*7.33e-3)));
% Zry: strain rate sensitivity exponent MATPRO
  clad.m = @(T) (T<=730).*0.02 + (T>730 & T<=900) .* (20.63172161 - 0.07704552983*T + 9.504843067e-05*T.^2 - 3.860960716e-08*T.^3) + (T>900) .* (-6.47e-02 + T * 2.203e-04);
% Zry: strain hardening exponent MATPRO
  clad.n = @(T) (T<1099.0772).*(-9.490e-2 + T.*(1.165e-3 + T.*(-1.992e-6 + T*9.588e-10))) + (T>=1099.0772 & T<1600).*(-0.22655119 + 2.5e-4*T) + (T>=1600).*0.17344880;

end