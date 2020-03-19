function status = writeResults(time,y,flag)

  global g fuel gap ingas clad cool fileID iTimeStep

  switch flag
      case 'init'
        fprintf('time =%12.5e(s)\n',time(1)); % print time
        if ~clad.fail
           fileID=fopen('results.m','w'); % open file for writing output
           fprintf(fileID,'function s=results\n\n'); % this file is a matlab function which returns the structure s
           iTimeStep = 0;
        else
           fileID=fopen('results.m','a+'); % open file for appending output
        end

      case 'done'
      % Stop stopwatch
        elapsedTime = toc(g.timerValue);
        fprintf(fileID,'%% elapsedTime = %12.5e s / %12.5e min / %12.5e hrs;\n\n',elapsedTime, elapsedTime/60, elapsedTime/3600);
        fclose(fileID); % close the output file

      otherwise
        if y(end) ~= y(end-1)
           fprintf('time =%12.5e(s) ',time(end));
           fprintf('cool.regime = '); % print regime
           fprintf('%1i',cool.regime);
           fprintf('\n');

           fprintf(fileID,'%s ',repmat('%',1,200));
           fprintf(fileID,'\n');
           iTimeStep = iTimeStep + 1;
           fprintf(fileID,'iTimeStep = %i;\n',iTimeStep);
           fprintf(fileID,'s.time(iTimeStep)=%12.5e; %%s\n',time(end));
           fprintf(fileID,'s.ingasP(iTimeStep)=%12.5e; %%MPa\n',ingas.p(end));

		   for i=1:fuel.nr
               fprintf(fileID,'s.fuel.T(:,%2i,iTimeStep)= [',i); 
               fprintf(fileID,'%9.2f ',fuel.T(:,i)-273.15); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% fuel temperature (C)\n');
		   end

           fprintf(fileID,'\n');

		   for i=1:clad.nr
               fprintf(fileID,'s.clad.T(:,%2i,iTimeStep)= [',i);
               fprintf(fileID,'%9.2f ',clad.T(:,i)-273.15); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% clad temperature (C)\n');
		   end

           fprintf(fileID,'\n');

           fprintf(fileID,'s.cool.T(:,iTimeStep)=    ['); 
           fprintf(fileID,'%9.2f ',[cool.pro.T]-273.15);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% mixture temperature (C)\n');

           fprintf(fileID,'s.cool.TS(:,iTimeStep)=   ['); 
           fprintf(fileID,'%9.2f ',[cool.pro.Tsat]-273.15);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% saturation temperature (C)\n');

           fprintf(fileID,'s.cool.TCHF(:,iTimeStep)= ['); 
           fprintf(fileID,'%9.2f ',cool.TCHF-273.15);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% critical heat flux temperature (C)\n');

           fprintf(fileID,'s.cool.regm(:,iTimeStep)= ['); 
           fprintf(fileID,'%9i ',cool.regime);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% flow regime (-)\n');

           fprintf(fileID,'\n');

           fprintf(fileID,'s.cool.p(:,iTimeStep)=    ['); 
           fprintf(fileID,'%9.6f ',cool.p);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% pressure (MPa)\n');

           fprintf(fileID,'s.cool.h(:,iTimeStep)=    [');
           fprintf(fileID,'%9.2f ',cool.h/1e3);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% enthalpy (kJ/kg)\n');

           fprintf(fileID,'s.cool.vel(:,iTimeStep)=  ['); 
           fprintf(fileID,'%9.5f ',cool.vel);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% mixture velocity (m/s)\n');

           fprintf(fileID,'s.cool.mdot(:,iTimeStep)= ['); 
           fprintf(fileID,'%9.5f ',cool.mdot);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% mixture flowrate (kg/s)\n');

           fprintf(fileID,'s.cool.genrG(:,iTimeStep)=['); 
           fprintf(fileID,'%9.5f ',cool.gamma);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% gas generation rate (kg/s)\n');

           fprintf(fileID,'s.cool.void(:,iTimeStep)= ['); 
           fprintf(fileID,'%9.5f ',[cool.pro.void]);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% void (-)\n');

           fprintf(fileID,'s.cool.x(:,iTimeStep)=    ['); 
           fprintf(fileID,'%9.5f ',[cool.pro.x]);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% equilibrium quality (-)\n');

           fprintf(fileID,'\n');

           fprintf(fileID,'s.fuel.r(:,iTimeStep)=    ['); 
           fprintf(fileID,'%9.5f ',fuel.r*1000);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% fuel outer radius (mm)\n');

           fprintf(fileID,'s.clad.rIn(:,iTimeStep)=  ['); 
           fprintf(fileID,'%9.5f ',clad.r(:,1)*1000);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% clad inner radius (mm)\n');

           fprintf(fileID,'s.clad.rOut(:,iTimeStep)= ['); 
           fprintf(fileID,'%9.5f ',clad.r(:,end)*1000);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% clad outer radius (mm)\n');

           fprintf(fileID,'\n');

           for i=1:clad.nr
               fprintf(fileID,'s.sigr(:,%2i,iTimeStep)=   [',i);
               fprintf(fileID,'%9.2f ',clad.sig{1}(:,i)); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% radial stress (MPa)\n');
		   end

           fprintf(fileID,'\n');

           for i=1:clad.nr
               fprintf(fileID,'s.sigh(:,%2i,iTimeStep)=   [',i);
               fprintf(fileID,'%9.2f ',clad.sig{2}(:,i)); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% hoop stress (MPa)\n');
		   end

           fprintf(fileID,'\n');

           for i=1:clad.nr
               fprintf(fileID,'s.sigz(:,%2i,iTimeStep)=   [',i);
               fprintf(fileID,'%9.2f ',clad.sig{3}(:,i)); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% axial stress (MPa)\n');
		   end

           fprintf(fileID,'\n');

           fprintf(fileID,'s.sigI(:,iTimeStep)=      ['); 
           fprintf(fileID,'%9.2f ',clad.sigI);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% clad engineering stress (MPa)\n');

           fprintf(fileID,'s.sigB(:,iTimeStep)=      ['); 
           fprintf(fileID,'%9.2f ',clad.sigB(clad.Tavg));
           fprintf(fileID,']; ');
           fprintf(fileID,'%% clad burst stress (MPa)\n');

           fprintf(fileID,'\n');

           fprintf(fileID,'s.gap.dr(:,iTimeStep)=    ['); 
           fprintf(fileID,'%9.2f ',gap.dr*1e6);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% gap width (mkm)\n');

           fprintf(fileID,'s.gap.h(:,iTimeStep)=     ['); 
           fprintf(fileID,'%9.2f ',gap.h);
           fprintf(fileID,']; ');
           fprintf(fileID,'%% gap conductance (W/m2K)\n');

           fprintf(fileID,'\n');

		   for i=1:clad.nr
               fprintf(fileID,'s.epsPeff(:,%2i,iTimeStep)=[',i);
               fprintf(fileID,'%9.2f ',clad.epsPeff(:,i)*100); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% effective plastic strain(%%)\n');
		   end

           fprintf(fileID,'\n');

		   for i=1:clad.nr
               fprintf(fileID,'s.epsrP(:,%2i,iTimeStep)=  [',i);
               fprintf(fileID,'%9.2f ',clad.epsP{1}(:,i)*100); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% radial plastic strain(%%)\n');
		   end

           fprintf(fileID,'\n');

		   for i=1:clad.nr
               fprintf(fileID,'s.epshP(:,%2i,iTimeStep)=  [',i);
               fprintf(fileID,'%9.2f ',clad.epsP{2}(:,i)*100); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% hoop plastic strain(%%)\n');
		   end

           fprintf(fileID,'\n');

		   for i=1:clad.nr
               fprintf(fileID,'s.epszP(:,%2i,iTimeStep)=  [',i);
               fprintf(fileID,'%9.2f ',clad.epsP{3}(:,i)*100); 
               fprintf(fileID,']; ');
               fprintf(fileID,'%% axial plastic strain(%%)\n');
		   end
        end
  end
  status = 0;
end