% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% First download the GENDF files for the required isotopes from the 
% open-access IAEA website:
% https://www-nds.iaea.org/ads/adsgendf.html and then processed them
% with the function convertGXStoCSV.
%
% The function scans the current folder for the files with the
% extension .GXS which are microscopic cross sections in 421 energy
% group structure in the GENDF format and convert them to the CSV
% ("Comma Separated Values") file format readable by Excel and
% therefore easily readable by Matlab. Note that semicolons are used
% in the script instead of commas, therefore, check that your Windows
% regional settings are set to use semicolons instead of commas as the
% list separator symbol.


  function convertGXStoCSV

% structure filesGXS contains information about all GXS files in the current directory
  filesGXS = dir('*.GXS');

% loop over all GXS files
  for i = 1:length({filesGXS.name})

    % define the name of the file without extension
      [~,nameOnly,~] = fileparts(filesGXS(i).name);

    % if the corresponding CSV file does not exist
      if exist([nameOnly,'.CSV'],'file') ~= 2

         fprintf('Convert %s to %s. Wait',[nameOnly,'.GXS'],[nameOnly,'.CSV']);
         fprintf(repmat('.',1,50));

       % open GENDF file for reading
         fdGXS = fopen([nameOnly,'.GXS']);
       % open CSV file for writing
         fdCSV = fopen([nameOnly,'.CSV'],'w');
       % counter of the bytes read
         bytesRead = 0;

       % cycle till the end of file
         while ~feof(fdGXS)

             % read the line
               str = fgets(fdGXS);
               bytesRead = bytesRead + length(str);
               if bytesRead/filesGXS(i).bytes > 0.02
                  bytesRead = 0;
                  fprintf ('\b \b');
               end

               ii = 1;
               for iii = 1:6
                   if str(ii+8) == '+' || str(ii+8) == '-' && str(ii+7) ~= ' '
                    % insert "E" (1.0+01 --> 1.0E+01)
                      str1 = [str(ii:ii+7),'E',str(ii+8:ii+10)];

                   elseif str(ii+9) == '+' || str(ii+9) =='-' && str(ii+8) ~= ' '
                    % insert "E" (1.0+01 --> 1.0E+01)
                      str1 = [str(ii:ii+8),'E',str(ii+9:ii+10)];

                   else
                      str1 = [' ',str(ii:ii+10)];

                   end

                 % write the line inserting semicolons
                   fprintf(fdCSV,' %s;',str1);
                   ii = ii + 11;
               end
                  
               fprintf(fdCSV,'%s;',str(67:70));
               fprintf(fdCSV,'%s;',str(71:72));
               fprintf(fdCSV,'%s;',str(73:75));
               fprintf(fdCSV,'%s\n',str(76:80));
               
         end
         fclose(fdGXS);
         fclose(fdCSV);
         fprintf(' Done.\n');
      end
  end