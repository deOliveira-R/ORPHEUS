% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function scans the current folder for the files with the extension
% .CSV which are microscopic cross sections in 421 energy group structure
% in the CSV ("Comma-Separated Value") format, obtained from the GENDF
% format, and convert them to the Matlab format.
%
% For every temperature available in the GENDF file, separate Matlab
% file is created.
%
% The Matlab format file for an isotope includes: 
% - atomic weight (amu);
% - number of energy groups (=421);
% - energy group boundaries (eV);
% - set of background cross sections, sigma-zeros (b)
% - temperature (K);
% - radiative capture cross sections (b) for each sigma-zero;
% - (n,alfa) reaction cross sections (b) for each sigma-zero;
% - (n,2n) matrix (b) for the first Legendre component;
% - scattering matrix (b) for the three Legendre components and each sigma-zero;
% - fission cross sections (b) for each sigma-zero;
% - nubar (-);
% - fission spectrum (-);
% - total cross sections (b) for each sigma-zero (b);
%
% Note that (n,3n), etc. reactions are NOT included in the current
% version.
%
% The GENDF files can be downloaded from the open-access IAEA website:
% https://www-nds.iaea.org/ads/adsgendf.html and processed with the function 
% convertGXStoCSV.


function convertCSVtoM

  filesCSV = dir('*.CSV'); % structure filesCSV contains information about all CSV files in the current directory
  for i = 1:length({filesCSV.name}) % loop over all CSV files
  
      [~,nameOnly,~] = fileparts(filesCSV(i).name); % find the name of the file without extension
      
      fprintf('Import data for %s and check if Matlab files for all temperatures are already available. ',nameOnly);
      m = importdata([nameOnly,'.CSV'],';'); % load CSV file into matrix m
      fprintf('Done.\n');
      nRow = size(m,1); % number of rows
         
    % Find number of temperatures and values of temperatures using mf=1 and mt=451
      nTemp = 0; % number of temperatures
      for iRow = 1:nRow
          if m(iRow,8) == 1 && m(iRow,9) == 451 && m(iRow,10) == 2
             nTemp = nTemp + 1; % number of temperatures
             temp(nTemp) = m(iRow,1); % vector of temperatures
          end
      end

      for iTemp = 1:nTemp % loop over temperatures

          if temp(iTemp) < 1000
             isoName = ['micro_',nameOnly,'__',num2str(round(temp(iTemp))),'K']; % name of the file with a temperature index
          else
             isoName = ['micro_',nameOnly,'_',num2str(round(temp(iTemp))),'K']; % name of the file with a temperature index
          end

          if exist([isoName,'.m'],'file') ~= 2 % if the corresponding Matlab file does not exist
             
             fdM = fopen([isoName,'.m'],'w'); % open the file for writing
           % Make a header for the file to be created with important parameters for
           % which the macroscopic cross sections were generated
             fprintf(fdM,'%% ---------------------------------------------------------\n');
             fprintf(fdM,'%% Matlab-based Open-source Reactor Physics Education System\n');
             fprintf(fdM,'%% ---------------------------------------------------------\n');
             fprintf(fdM,'%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.\n');
             fprintf(fdM,'%%\n');
             fprintf(fdM,'%% Microscopic cross sections for %s at %s K in 421 energy group structure\n', nameOnly,num2str(round(temp(iTemp))));
             fprintf(fdM,'%% generated from GENDF files downloaded from the IAEA website https://www-nds.iaea.org/ads/adsgendf.html\n\n');
             fprintf(fdM,'function s = %s\n\n',isoName); % write the first line of the function
             
             fprintf(fdM,'fprintf(''Read microscopic cross sections of %s.\\n'');\n\n',isoName);
             
             fprintf(fdM,'%% atomic weight (amu)\n');
             fprintf(fdM,'s.aw = %13.6e;\n\n', m(2,2)*1.008664916);
             
             ng = 421;
             fprintf(fdM,'%% number of energy groups\n');
             fprintf(fdM,'ng = %3i;\n\n', ng);
             
             nSig0 = m(2,4);
             a = extractNwords(1+nSig0+(ng+1), 4, m);
             fprintf(fdM,'%% energy group boundaries (eV)\n');
             fprintf(fdM,'s.eg = ['); 
             fprintf(fdM,'%12.5e ', a(2+nSig0 : 2+nSig0+ng)); 
             fprintf(fdM,'];\n\n');

             fprintf(fdM,'%% sigma-zeros (b)\n');
             fprintf(fdM,'s.sig0 = ['); 
             fprintf(fdM,'%12.5e ', a(2 :2+nSig0-1)); 
             fprintf(fdM,'];\n');
             fprintf(fdM,'nSig0 = %2i;\n\n', nSig0); 

             fprintf(fdM,'%% temperature (K)\n');
             fprintf(fdM,'s.temp = %9.2f;\n\n', temp(iTemp));

% (n,gamma)
             fprintf('Convert %s.CSV to %s.m: mf=3 mt=102 radiative capture\n',nameOnly,isoName);
             sigC = extract_mf3(102, iTemp, m); % Extract mf=3 mt=102 (radiative capture cross sections)
             nSig0C = size(sigC,1);
             fprintf(fdM,'%% radiative capture cross section (b) for %2i sigma-zero(s)\n',nSig0C);
             for iSig0 = 1:nSig0C
                 fprintf(fdM,'s.sigC(%2i,:)=[',iSig0);
                 fprintf(fdM,'%13.6e ',sigC(iSig0,1:ng));
                 fprintf(fdM,'];\n');
             end
             if nSig0C == 1 && nSig0 > 1
                fprintf(fdM,'sigC(2:nSig0,:) = sigC(1,:);\n');
             end
             fprintf(fdM,'\n');

% (n,alfa)
             fprintf('Convert %s.CSV to %s.m: mf=3 mt=107 (n,alfa)\n',nameOnly,isoName);
             sigL = extract_mf3(107, iTemp, m); % Extract mf=3 mt=107 (production of an alfa particle)
             if sigL == 0
                sigL = zeros(nSig0,ng);
                fprintf(fdM,'%% (n,alfa) cross section (b)\n');
                fprintf(fdM,'s.sigL = zeros(nSig0,ng);\n\n');
             else
                nSig0L = size(sigL,1);
                fprintf(fdM,'%% (n,alfa) cross sections (b) for %2i sigma-zero(s)\n',nSig0L);
                for iSig0 = 1:nSig0L
                    fprintf(fdM,'s.sigL(%2i,:)=[',iSig0);
                    fprintf(fdM,'%13.6e ',sigL(iSig0,1:ng));
                    fprintf(fdM,'];\n');
                end
                if nSig0L == 1 && nSig0 > 1
                   fprintf(fdM,'s.sigL = repmat(s.sigL,nSig0,1);\n');
                   sigL = repmat(sigL,nSig0,1);
                end
                fprintf(fdM,'\n');
             end

% (n,2n)
             [ifrom2, ito2, sig2] = extract_mf6(16, iTemp, m); % Extract mf=6 mt=16 ((n,2n) matrix)
             fprintf(fdM,'%% (n,2n) matrix for 1 Legendre components\n');
			 if ifrom2(1) == 0
                isn2n = 0;
                fprintf(fdM,'s.sig2=zeros(ng,ng);\n');
                fprintf(fdM,'\n');
             else
                isn2n = 1;
                fprintf('Convert %s.CSV to %s.m: mf=6 mt=16 (n,2n) reaction\n',nameOnly,isoName);
                fprintf(fdM,'ifrom=[');fprintf(fdM,blanks(18));fprintf(fdM,'%13i ',ifrom2);fprintf(fdM,'];\n');
                fprintf(fdM,'ito  =[');fprintf(fdM,blanks(18));fprintf(fdM,'%13i ',ito2);fprintf(fdM,'];\n');
                fprintf(fdM,'s.sig2=sparse(ifrom,ito,[');
                fprintf(fdM,'%13.6e ',sig2{1+0,1});
                fprintf(fdM,'],ng,ng);\n');
                fprintf(fdM,'\n');
			 end

% (n,n')
             igThresh = 95; % last group of thermal energy (e = 4 eV)

             fprintf('Convert %s.CSV to %s.m: mf=6 mt=2 elastic scattering\n',nameOnly,isoName);
             [ifromE, itoE, sigE] = extract_mf6(2, iTemp, m); % Extract mf=6 mt=2 (elastic scattering matrix)
             nLgn = size(sigE,1)-1;
             for jLgn = 0:nLgn
                 for iSig0 = 1:nSig0
                     for ii = 1:length(ifromE)
                         if ifromE(ii) <= igThresh
                            sigE{1+jLgn,iSig0}(ii) = 0;
                         end
                     end
                     sigS{1+jLgn,iSig0} = sparse(ifromE, itoE, sigE{1+jLgn,iSig0}+1e-30, ng, ng);
                 end
             end

             for ii=51:91
                 [ifromI, itoI, sigI] = extract_mf6(ii, iTemp, m); % Extract mf=6 mt=51 ... 91 (inelastic scattering matrix)
                 if ifromI(1) > 0
                    fprintf('Convert %s.CSV to %s.m: mf=6 mt=%2i inelastic scattering\n',nameOnly,isoName,ii);
                    nLgn = size(sigI,1)-1;
                    for jLgn = 0:nLgn
                        for iSig0 = 1:nSig0
                            sigS{1+jLgn,iSig0} = sigS{1+jLgn,iSig0} + sparse(ifromI, itoI, sigI{1+jLgn,1}+1e-30, ng, ng);
                        end
                    end
                 end
             end

             if isoName(1:11) == 'micro_H_001' %#ok<*STCMP>
                fprintf('Convert %s.CSV to %s.m: mf=6 mt=222 thermal scattering for hydrogen binded in water\n',nameOnly,isoName);
                [ifromI, itoI, sigI] = extract_mf6(222, iTemp, m); % Extract mf=6 mt=222 thermal scattering for hydrogen binded in water
                nLgn = size(sigI,1)-1;
                for jLgn = 0:nLgn
                    for iSig0 = 1:nSig0
                        sigS{1+jLgn,iSig0} = sigS{1+jLgn,iSig0} + sparse(ifromI, itoI, sigI{1+jLgn,1}+1e-30, ng, ng);
                    end
                end
             else
                fprintf('Convert %s.CSV to %s.m: mf=6 mt=221 free gas thermal scattering\n',nameOnly,isoName);
                [ifromI, itoI, sigI] = extract_mf6(221, iTemp, m); % Extract mf=6 mt=221 free gas thermal scattering
                nLgn = size(sigI,1)-1;
                for jLgn = 0:nLgn
                    for iSig0 = 1:nSig0
                        sigS{1+jLgn,iSig0} = sigS{1+jLgn,iSig0} + sparse(ifromI, itoI, sigI{1+jLgn,1}+1e-30, ng, ng);
                    end
                end
             end

             fprintf(fdM,'%% scattering matrix for 3 Legendre components for %2i sigma-zeros\n',nSig0);
             for jLgn = 0:2
                 notYetPrinted = true;
                 for iSig0 = 1:nSig0
                     [ifrom,ito,sigS_] = find(sigS{1+jLgn,iSig0}); % Find indices and values of nonzero elements
                   % ifrom and ito are printed only once
                     if notYetPrinted
                        fprintf(fdM,'ifrom=[');fprintf(fdM,blanks(26));fprintf(fdM,'%13i ',ifrom);fprintf(fdM,'];\n');
                        fprintf(fdM,'ito  =[');fprintf(fdM,blanks(26));fprintf(fdM,'%13i ',ito);fprintf(fdM,'];\n');
                        notYetPrinted = false;
                     end
                     fprintf(fdM,'s.sigS{1+%1i,%2i}=sparse(ifrom,ito,[',jLgn,iSig0); 
                     fprintf(fdM,'%13.6e ',sigS_);
                     fprintf(fdM,'],ng,ng);\n');
                 end
                 fprintf(fdM,'\n');
             end

% (n,fis)             
             sigF = extract_mf3(18, iTemp, m); % Extract mf=3 mt=18 (fission cross sections)
             if size(sigF) == 1
                sigF = zeros(nSig0,ng);                
                fprintf(fdM,'s.fissile = 0;\n\n');
                
                fprintf(fdM,'%% fission cross sections (b)\n');
                fprintf(fdM,'s.sigF = zeros(nSig0,ng);\n\n');
                
                fprintf(fdM,'%% nubar\n');
                fprintf(fdM,'s.nubar = zeros(nSig0,ng);\n\n');
                
                fprintf(fdM,'%% fission spectrum\n');
                fprintf(fdM,'s.chi = zeros(nSig0,ng);\n\n');
             else   
                fprintf('Convert %s.CSV to %s.m: mf=3 mt=18 fission\n',nameOnly,isoName);
                fprintf(fdM,'s.fissile=1;\n\n');
                
                nSig0F = size(sigF,1);
                fprintf(fdM,'%% fission cross sections (b) for %2i sigma-zero(s)\n',nSig0F);
                for iSig0 = 1:nSig0F;
                    fprintf(fdM,'s.sigF(%2i,:)=[',iSig0);
                    fprintf(fdM,'%13.6e ',sigF(iSig0,1:ng));
                    fprintf(fdM,'];\n');
                end
                fprintf(fdM,'\n');

                nubar = extract_mf3(452, iTemp, m); % Extract mf=3 mt=452 (total nubar)
                fprintf('Convert %s.CSV to %s.m: mf=3 mt=452 total nubar\n',nameOnly,isoName);
                nSig0nu = size(nubar,1);
                fprintf(fdM,'%% total nubar for %2i sigma-zero(s)\n',nSig0nu);
                for iSig0 = 1:nSig0nu
                    fprintf(fdM,'s.nubar(%2i,:)=[',iSig0);
                    fprintf(fdM,'%13.6e ',nubar(iSig0,1:ng));
                    fprintf(fdM,'];\n');
                end
                fprintf(fdM,'\n');

% chi
                fprintf('Convert %s.CSV to %s.m: mf=6 mt=18 fission spectrum\n',nameOnly,isoName);
                iRow=1;
                while m(iRow,8) ~= 6 || m(iRow,9) ~= 18 % find fission spectrum
                   iRow = iRow + 1;     
                end
                iRow = iRow + 1;
                ig2lo = m(iRow,4); % index to lowest nonzero group
                nw = m(iRow,5); % number of words to be read
                iRow = iRow + 1;     
                [a,~] = extractNwords(nw, iRow, m); % read nw words in vector a
                for iii = 1:ig2lo-1
                    chi(iii)=0.0;
                end
                for iii = 1:nw
                    chi(iii+ig2lo-1) = a(iii);
                end
                
                fprintf(fdM,'%% fission spectrum\n');
                fprintf(fdM,'s.chi=[');
                fprintf(fdM,'%13.6e ',chi(1:ng)/sum(chi));
                fprintf(fdM,'];\n\n');
             end
             
% total
           % Calculate total cross sections (note that mf=3 mt=1 does not
           % include upscatters).
             fprintf(fdM,'%% total cross sections (b) for %2i sigma-zeros\n',nSig0);
             for iSig0 = 1:nSig0
                 sigT(iSig0,:) = sigC(iSig0,:) + sigF(iSig0,:) + sigL(iSig0,:) + sum(sigS{1+0,iSig0},2)';
                 if isn2n
                    sigT(iSig0,:) = sigT(iSig0,:) + sum(sig2{1+0,1},2)';
                 end
                 fprintf(fdM,'s.sigT(%2i,:)=[',iSig0);
                 fprintf(fdM,'%13.6e ',sigT(iSig0,:));
                 fprintf(fdM,'];\n');
             end
             fprintf(fdM,'\n');
          end % of 'if the file does not exist'
      end % of loop over temperatures
  end
end
  
  
% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function searches matrix m for cross sections sig from file mf=3 for
% reaction mt and temperature ntt and and returns sig(ng,nSig0), where ng
% is the number of energy groups and nSig0 is the the number of
% sigma-zeros.
%
%
function sig = extract_mf3(mt, ntt, m)

  nRow = size(m,1); % number of rows
  nTemp = 0; % number of temperatures
  iRowFound = 0;
  for iRow = 1:nRow
      if m(iRow,8) == 3 && m(iRow,9) == mt && m(iRow,10) == 1 % find the row with mf=3 and required mt
         nTemp = nTemp + 1; % number of temperatures
         if nTemp == ntt
            iRowFound = iRow+1;
            break;
         end
      end
  end

  if iRowFound > 0 % there is mf=3 and required mt for this isotope
     nSig0 = m(iRowFound-1,4); % number of sigma-zeros
     nLgn = m(iRowFound-1,3); % number of Legendre components
     iRow = iRowFound + 1;
     while m(iRow,8) == 3 && m(iRow,9) == mt
         ig = m(iRow-1,6);
         [a, iRowNew] = extractNwords(nSig0*nLgn*2, iRow, m); 
         sig(1:nSig0,ig) = a(nSig0*nLgn+(1:nSig0));
         iRow = iRowNew + 2;
     end
  else
     sig = 0;
  end
end






% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads cross sections from file 6 for reaction mt and 
% temperature index ntt from matrix m and returns the 2D cell matrix 
% sig{nLgn,nSig0}(nonz) with two vectors ifrom(nonz) and ito(nonz), where
% nLgn is the number of Legendre components, nSig0 is the number of
% sigma-zeros and nonz is the number of nonzeros.
 
 function [ifrom,ito,sig] = extract_mf6(mt, ntt, m)
 
  iRow = 1; % row number
  nTemp = 0; % number of temperatures
  ifrom = 0; % index of group 'from'
  ito = 0; % index of group 'to'
  
  while m(iRow,7) ~= -1 % up to the end
      if m(iRow,8) == 6 && m(iRow,9) == mt % find the row with mf=6 & mt
         if m(iRow,10) == 1  % this is the first line of mf=6 & mt: initialize
            nonz = 0; % number of nonzeros
            nLgn = m(iRow,3); % number of Legendre components
            nSig0 = m(iRow,4); % number of sigma-zeros

            iRow = iRow + 1;
            nTemp = nTemp + 1; % temperature index
         end
         ng2 = m(iRow,3); % number of secondary positions
         ig2lo = m(iRow,4); % index to lowest nonzero group
         nw = m(iRow,5); % number of words to be read
         ig = m(iRow,6); % current group index

         iRow = iRow + 1;     
         [a,iRowNew] = extractNwords(nw,iRow,m); % extract nw words in vector a
         iRow = iRowNew;

         if nTemp == ntt
            k = nLgn*nSig0; % the first nLgn*nSig0 words are flux -- skip.
            for iTo = ig2lo : ig2lo+ng2-2
                nonz = nonz + 1;
                ifrom(nonz) = ig;
                ito(nonz) = iTo;
                for iSig0 = 1:nSig0
                    for iLgn = 1:nLgn
                        k = k + 1;
                        sig{iLgn,iSig0}(nonz) = a(k);
                    end
                    if nLgn == 1
                       sig{1+1,iSig0}(nonz) = 0;
                       sig{1+2,iSig0}(nonz) = 0;
                    end
                end
            end    
         end
      end
      iRow = iRow + 1;
  end
  if nTemp == 0
      sig = 0;
  end
end



% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads n words from row iRow of matrix m and returns them in
% vector a together with the new row number iRowNew, i.e. the row where the
% last word was read.
%
%
 function [a, iRowNew] = extractNwords(n, iRow, m)
  
  k = 1; 
  iRowNew = iRow;
  for ii = 1:fix(n/6) %read lines with 6 words each
      for jj = 1:6
          a(k) = m(iRowNew,jj);
          k = k + 1;
      end
      iRowNew = iRowNew + 1;
  end
  
  if (n - fix(n/6)*6) == 0
     iRowNew = iRowNew - 1;
  end
  
  for jj = 1:(n-fix(n/6)*6) %read the last line with less than 6 words
      a(k) = m(iRowNew,jj); 
      k = k + 1;
  end
 end
