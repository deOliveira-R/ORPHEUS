  fprintf('Plot the results...\n');

%--------------------------------------------------------------------------
% Plot the results
  s = results;

  f = figure('visible','off');
  plot(s.time(:),s.pow,'-b');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Relative power (-)');
  saveas(f, 'Fig_pow.pdf');

  f = figure('visible','off');
  plot(s.time(:),s.reacDop,'-r', ...
	   s.time(:),s.reacCool,'-b', ...
	   s.time(:),s.reac,'-k');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Reactivity (pcm)');
  legend('Doppler', 'coolant temperature', 'total', ...
         'Location','best');
  saveas(f, 'Fig_reac.pdf');

  f = figure('visible','off');
  plot(s.time(:),s.cool.void(end,:),'-b');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Void fraction (-)');
  saveas(f, 'Fig_void.pdf');

  f = figure('visible','off');
  plot(s.time(:),squeeze(s.fuel.T(end,1,:)),'-r', ...
       s.time(:),squeeze(s.clad.T(end,end,:)),'-g', ...
       s.time(:),s.cool.T(end,:),'-b', ...
       s.time(:),s.cool.TS(end,:),'--b',...
       s.time(:),s.cool.TCHF(end,:),'--g');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Temperature at outlet (C)');
  legend('centre of fuel','outer surface of clad', 'coolant', 'coolant saturation', 'CHF',...
         'Location','best');
  saveas(f, 'Fig_T.pdf');

  f = figure('visible','off');
  plot(s.time(:),squeeze(s.epsPeff(end,1,:)),'-r', ...
       s.time(:),squeeze(s.epsrP(end,1,:)),'-g', ...
       s.time(:),squeeze(s.epshP(end,1,:)),'-b', ...
       s.time(:),squeeze(s.epszP(end,1,:)),'-k');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Plastic strain at inner clad surface (%)');
  legend('effective', 'radial', 'hoop', 'axial',...
         'Location','best');
  saveas(f, 'Fig_epsP.pdf');

  f = figure('visible','off');
  plot(s.time,s.ingasP,'-r', ...
       s.time,s.cool.p(end,:),'-b');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Pressure (MPa)');
  legend('inner gas', 'coolant',...
         'Location','best');
  saveas(f, 'Fig_pressure.pdf');

  f = figure('visible','off');
  plot(s.time,s.fuel.r(end,:),'-r', ...
       s.time,s.clad.rIn(end,:),'-b', ...
       s.time,s.clad.rOut(end,:),'-g');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Radius (mm)');
  legend('outer fuel', 'inner clad', 'outer clad', ...
         'Location','best');
  saveas(f, 'Fig_radii.pdf');

  f = figure('visible','off');
  plot(s.time,squeeze(s.sigr(2,1,:)),'-r', ...
       s.time,squeeze(s.sigh(2,1,:)),'-g', ...
       s.time,squeeze(s.sigz(2,1,:)),'-b');
  ylim(ylim);
  grid on;
  xlabel('Time (s)');
  ylabel('Clad stress (MPs)');
  legend('Engineering stress', 'Burst stress', ...
         'Location','best');
  saveas(f, 'Fig_sig.pdf');
