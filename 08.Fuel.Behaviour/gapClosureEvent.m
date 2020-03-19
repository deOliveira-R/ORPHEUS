function [value,isterminal,direction] = gapClosureEvent(t,y)

% Global structure with geometry and thermal-hydraulic parameters 
  global gap

  value = gap.dr - gap.rough;
  isterminal = 1;
  direction = 0;

end