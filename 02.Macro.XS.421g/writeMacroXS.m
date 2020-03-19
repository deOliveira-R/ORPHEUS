% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function writes all group macroscopic cross sections from structure s
% to file matName.m

function writeMacroXS(s, matName)


  fprintf('Write macroscopic cross sections to the file: %s',[matName,'.m'])
  
  fdM=fopen([matName,'.m'],'w');
  
  fprintf(fdM,'function s = %s\n\n',matName);

  for i = 1:length(s.header)
      fprintf(fdM,'s.header{%2i}=''%s'';\n',i,s.header{i});
  end
  fprintf(fdM,'\n\n');
  
  fprintf(fdM,'%% atomic weight (amu)\n');
  fprintf(fdM,'s.aw = %12.5e;\n\n',s.aw);
  
  fprintf(fdM,'%% density of material (g/cm3)\n');
  fprintf(fdM,'s.den = %12.5e;\n\n',s.den);

  fprintf(fdM,'%% temperature (K)\n');
  fprintf(fdM,'s.temp = %9.2f;\n\n', s.temp);

  fprintf(fdM,'%% number of energy groups\n');
  fprintf(fdM,'s.ng = %3i;\n\n',s.ng);
  
  fprintf(fdM,'%% energy group boundaries (eV)\n');
  fprintf(fdM,'s.eg = [');
  fprintf(fdM,'%12.5e ',s.eg);
  fprintf(fdM,'];\n\n');

  fprintf(fdM,'%% isotope names\n');
  fprintf(fdM,'s.isoName = { ...\n');
  for i = 1:length(s.isoName)
      fprintf(fdM,'''%s'' ...\n',s.isoName{i});
  end
  fprintf(fdM,'};\n\n');
  
  fprintf(fdM,'%% isotope number densities (1/cm3)\n');
  fprintf(fdM,'s.numDen = [ ...\n');
  for i = 1:length(s.numDen)
      fprintf(fdM,'%12.5e ...\n',s.numDen(i));
  end
  fprintf(fdM,'];\n\n');
  
  if isfield(s,'sig0')
     fprintf(fdM,'%% background cross sections, sigma-zeros, for %2i isotopes (barn)\n',length(s.isoName));
     for i = 1:size(s.sig0,1)
         fprintf(fdM,'s.sig0(%2i,:) = [',i);
         fprintf(fdM,'%13.6e ',s.sig0(i,1:s.ng));
         fprintf(fdM,'];\n');
     end
     fprintf(fdM,'\n');
  end
  
  fprintf(fdM,'%% radiative capture cross section(1/cm)\n');
  fprintf(fdM,'s.SigC=[');
  fprintf(fdM,'%13.6e ',s.SigC(1:s.ng));
  fprintf(fdM,'];\n\n');
  
  fprintf(fdM,'%% (n,alfa) reaction cross section(1/cm)\n');
  fprintf(fdM,'s.SigL=[');
  fprintf(fdM,'%13.6e ',s.SigL(1:s.ng));
  fprintf(fdM,'];\n\n');

  fprintf(fdM,'%% scattering matrix for %1i Legendre components\n',size(s.SigS,2));
  for j=1:size(s.SigS,2)
      [ifrom,ito,Sig{j}] = find(s.SigS{j});
      fprintf(fdM,'ifrom=[');fprintf(fdM,'%13i ',ifrom);fprintf(fdM,'];\n');
      fprintf(fdM,'ito  =[');fprintf(fdM,'%13i ',ito);fprintf(fdM,'];\n');
      fprintf(fdM,'Sig  =[');fprintf(fdM,'%13.6e ',Sig{j});fprintf(fdM,'];\n');
      fprintf(fdM,'s.SigS{1+%1i} =sparse(ifrom,ito,Sig,s.ng,s.ng);\n\n',j-1);
  end
  
  fprintf(fdM,'%% (n,2n) matrix for 1 Legendre component\n');
  [ifrom,ito,Sig] = find(s.Sig2);
  fprintf(fdM,'ifrom=[');fprintf(fdM,'%13i ',ifrom);fprintf(fdM,'];\n');
  fprintf(fdM,'ito  =[');fprintf(fdM,'%13i ',ito);fprintf(fdM,'];\n');
  fprintf(fdM,'Sig  =[');fprintf(fdM,'%13.6e ',Sig);fprintf(fdM,'];\n');
  fprintf(fdM,'s.Sig2=sparse(ifrom,ito,Sig,s.ng,s.ng);\n\n');
  
  if s.SigP(1) ~= 0
     fprintf(fdM,'s.fissile=1;\n\n');
	 
     fprintf(fdM,'%% fission cross section(1/cm)\n');
     fprintf(fdM,'s.SigF=[');
     fprintf(fdM,'%13.6e ',s.SigF(1:s.ng));
     fprintf(fdM,'];\n\n');
	 
     fprintf(fdM,'%% production cross section(1/cm)\n');
     fprintf(fdM,'s.SigP=[');
     fprintf(fdM,'%13.6e ',s.SigP(1:s.ng));
     fprintf(fdM,'];\n\n');

     fprintf(fdM,'%% fission spectrum\n');
     fprintf(fdM,'s.chi=[');
     fprintf(fdM,'%13.6e ',s.chi(1:s.ng));
     fprintf(fdM,'];\n\n');
 
  else % s.SigP(1) == 0
     fprintf(fdM,'s.fissile=0;\n\n');
	 
     fprintf(fdM,'%% fission cross section(1/cm)\n');
     fprintf(fdM,'s.SigF=zeros(1,s.ng);\n\n');

     fprintf(fdM,'%% production cross section(1/cm)\n');
     fprintf(fdM,'s.SigP=zeros(1,s.ng);\n\n');

     fprintf(fdM,'%% fission spectrum\n');
     fprintf(fdM,'s.chi=zeros(1,s.ng);\n\n');
	 
  end
  
  fprintf(fdM,'%% total cross section(1/cm)\n');
  fprintf(fdM,'s.SigT=[');
  fprintf(fdM,'%13.6e ',s.SigT(1:s.ng));
  fprintf(fdM,'];\n\n');
  
  fprintf(fdM,'end');
  fclose(fdM);
  
  fprintf('. Done.\n')
end
