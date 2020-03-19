% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function fly calculates the angular group flux fiB in node B
% given the angular group flux fiA in node A, the source (qAB) of
% neutrons born between A and B and able to reach B (assumed constant
% between the nodes), the transport-corrected total cross section
% sigAB between the nodes and the distance sAB between the nodes.
%
% A      qAB, sigAB      B
% o----------------------o
% |<---------sAB-------->|

function fiB = fly(fiA, sigAB, sAB, qAB)

   fiB = fiA.*exp(-sigAB*sAB) + qAB.*(1 - exp(-sigAB*sAB))./sigAB;
   
end
