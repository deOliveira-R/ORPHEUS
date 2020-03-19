function status = writeResults(time,y,flag)

  global fdres iTimeStep fuel gap ingas clad
  
  switch flag
      case 'init'
        fprintf('time =%12.5e(s)\n',time(1)); % print time
        if gap.open
           fdres=fopen('results.m','w'); % open file for writing output
           fprintf(fdres,'function s=results\n\n'); % this file is a matlab function which returns the structure s
           iTimeStep = 0;
        else
           fdres=fopen('results.m','a+'); % open file for appending output
        end
        
      case 'done'
        if gap.clsd
           fprintf(fdres,'end');
        end
        fclose(fdres); % close the output file 
        
      otherwise
        if y(end) ~= y(end-1)
           fprintf('time =%10.3e(s) = %10.3e(min) = %10.3e(h) = %10.3e(d) = %12.5e(y); ',time(end),time(end)/60,time(end)/3600,time(end)/86400,time(end)/31536000);
           if gap.open
              fprintf('open gap width = %3.0f(mkm)\n', gap.dr*1e6);
           else
              fprintf('closed gap width = %3.0f(mkm)\n', gap.dr*1e6);
           end
           
           fprintf(fdres,'%s ',repmat('%',1,200));
           fprintf(fdres,'\n');

           iTimeStep = iTimeStep + 1;
           fprintf(fdres,'iTimeStep = %i;\n',iTimeStep);
           fprintf(fdres,'s.timeS(iTimeStep)=%12.5e; %%s\n',time(end));
           fprintf(fdres,'s.timeY(iTimeStep)=%12.5e; %%y\n',time(end)/31536000);
           fprintf(fdres,'s.LHGR(iTimeStep)=%12.5e; %%linear heat generation rate (W/m)\n',fuel.qLHGR);
           fprintf(fdres,'s.fuel.Bu =%12.5e; %%fuel burnup (MW-d/kgFuel)\n',fuel.Bu);
           fprintf(fdres,'s.fuel.innR(iTimeStep)=%12.5e; %%fuel inner radius (mm)\n',fuel.r(1)*1e3);
           fprintf(fdres,'s.fuel.outR(iTimeStep)=%12.5e; %%fuel outer radius (mm)\n',fuel.r(fuel.nr)*1e3);
           fprintf(fdres,'s.fuel.dz(iTimeStep)=%12.5e; %%fuel height (mm)\n',fuel.dz(1)*1e3);
           fprintf(fdres,'s.clad.innR(iTimeStep)=%12.5e; %%clad inner radius (mm)\n',clad.r(1)*1e3);
           fprintf(fdres,'s.clad.outR(iTimeStep)=%12.5e; %%clad outer radius (mm)\n',clad.r(clad.nr)*1e3);
           fprintf(fdres,'s.clad.dz(iTimeStep)=%12.5e; %%clad height (mm)\n',clad.dz(1)*1e3);
           fprintf(fdres,'s.gap.dr(iTimeStep)=%12.5e; %%gap width (mkm)\n',gap.dr*1e6);
           fprintf(fdres,'s.gap.open(iTimeStep)=%1i; %%gap open flag\n',gap.open);
           fprintf(fdres,'s.gap.clsd(iTimeStep)=%1i; %%gap closed flag\n',gap.clsd);
           fprintf(fdres,'s.gap.kGasMix(iTimeStep)=%12.5e; %% gap gas mixtire conductivity (W/m-K)\n',gap.kGasMix);
           fprintf(fdres,'s.gap.h(iTimeStep)=%12.5e; %% fuel-clad gap conductance (W/m2-K)\n',gap.h);
           fprintf(fdres,'s.gap.pcontact(iTimeStep)=%12.5e; %% fuel-clad contact pressure (MPa)\n',gap.pContact);
           fprintf(fdres,'s.ingas.p(iTimeStep)=%12.5e; %%inner gas pressure (MPa)\n',ingas.p);

           fprintf(fdres,'s.fuel.r(1:%2i,iTimeStep)=      [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.r*1e3); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel node radii (mm)\n');

           fprintf(fdres,'s.fuel.T(1:%2i,iTimeStep)=      [',fuel.nr); 
           fprintf(fdres,'%9.2f ',fuel.T-273.15); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel temperature (C)\n');

           fprintf(fdres,'s.fuel.sigVM(1:%2i,iTimeStep)=  [',fuel.nr); 
           fprintf(fdres,'%9.2f ',fuel.sigVM); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel Von Misses stress (MPa)\n');

           fprintf(fdres,'s.fuel.sigr(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.2f ',fuel.sig{1}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel radial stress (MPa)\n');

           fprintf(fdres,'s.fuel.sigh(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.2f ',fuel.sig{2}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel hoop stress (MPa)\n');

           fprintf(fdres,'s.fuel.sigz(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.2f ',fuel.sig{3}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel axial stress (MPa)\n');

           fprintf(fdres,'s.fuel.epsr(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.eps{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel radial strain (%%)\n');

           fprintf(fdres,'s.fuel.epsh(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.eps{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel hoop strain (%%)\n');

           fprintf(fdres,'s.fuel.epsz(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.eps{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel axial strain (%%)\n');

           fprintf(fdres,'s.fuel.epsrE(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsE{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel radial elastic strain (%%)\n');

           fprintf(fdres,'s.fuel.epshE(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsE{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel hoop elastic strain (%%)\n');

           fprintf(fdres,'s.fuel.epszE(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsE{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel axial elastic strain (%%)\n');

           fprintf(fdres,'s.fuel.epsT(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsT*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel linear thermal expansion (%%)\n');

           fprintf(fdres,'s.fuel.epsrC(1:%2i,iTimeStep)=  [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsC{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel creep radial strain (%%)\n');

           fprintf(fdres,'s.fuel.epshC(1:%2i,iTimeStep)=  [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsC{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel creep hoop strain (%%)\n');

           fprintf(fdres,'s.fuel.epszC(1:%2i,iTimeStep)=  [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsC{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel creep axial strain (%%)\n');

           fprintf(fdres,'s.fuel.epsS(1:%2i,iTimeStep)=   [',fuel.nr); 
           fprintf(fdres,'%9.5f ',fuel.epsS*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% fuel volumetric swelling deformation (%%)\n');

           fprintf(fdres,'s.clad.r(1:%2i,iTimeStep)=      [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.r*1e3); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad node radii (mm)\n');

           fprintf(fdres,'s.clad.T(1:%2i,iTimeStep)=      [',clad.nr); 
           fprintf(fdres,'%9.2f ',clad.T-273.15); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad temperature (C)\n');

           fprintf(fdres,'s.clad.sigVM(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.2f ',clad.sigVM); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad Von Misses stress (MPa)\n');

           fprintf(fdres,'s.clad.sigr(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.2f ',clad.sig{1}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad radial stress (MPa)\n');

           fprintf(fdres,'s.clad.sigh(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.2f ',clad.sig{2}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad hoop stress (MPa)\n');

           fprintf(fdres,'s.clad.sigz(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.2f ',clad.sig{3}); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad axial stress (MPa)\n');

           fprintf(fdres,'s.clad.epsr(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.eps{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad radial strain (%%)\n');

           fprintf(fdres,'s.clad.epsh(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.eps{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad hoop strain (%%)\n');

           fprintf(fdres,'s.clad.epsz(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.eps{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad axial strain (%%)\n');

           fprintf(fdres,'s.clad.epsrE(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsE{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad radial elastic strain (%%)\n');

           fprintf(fdres,'s.clad.epshE(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsE{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad hoop elastic strain (%%)\n');

           fprintf(fdres,'s.clad.epszE(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsE{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad axial elastic strain (%%)\n');

           fprintf(fdres,'s.clad.epsT(1:%2i,iTimeStep)=   [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsT*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad linear thermal expansion (%%)\n');

           fprintf(fdres,'s.clad.epsrC(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsC{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad creep radial strain (%%)\n');

           fprintf(fdres,'s.clad.epshC(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsC{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad creep hoop strain (%%)\n');

           fprintf(fdres,'s.clad.epszC(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsC{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad creep axial strain (%%)\n');

           fprintf(fdres,'s.clad.epsrP(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsP{1}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad plastic radial strain (%%)\n');

           fprintf(fdres,'s.clad.epshP(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsP{2}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad plastic hoop strain (%%)\n');

           fprintf(fdres,'s.clad.epszP(1:%2i,iTimeStep)=  [',clad.nr); 
           fprintf(fdres,'%9.5f ',clad.epsP{3}*100); 
           fprintf(fdres,']; ');
           fprintf(fdres,'%% clad plastic axial strain (%%)\n');
        end
  end
  status = 0;
