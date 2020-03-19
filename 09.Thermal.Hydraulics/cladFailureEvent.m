function [value,isterminal,direction] = cladFailureEvent(t,y)

% global structures
  global clad

  value = clad.sigB(clad.Tavg(end)) - clad.sigI(end);
  isterminal = 1;
  direction = 0;

end