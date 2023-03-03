function rMOON = moonpos(JDut1)
tUT1 = (JDut1-2451545.0)/36525.0;
lEC = 218.32 + 481267.8813*tUT1 + 6.29*sind(134.9+477198.85*tUT1) - 1.27*sind(259.2 - 413335.38*tUT1) + 0.66*sind(235.7+890534.23*tUT1) + 0.21*sind(269.9 + 954397.70*tUT1) - 0.19*sind(357.5+35999.05*tUT1) - 0.11*sind(186.6 + 966404.05*tUT1);

phiEC = 5.13*sind(93.3 + 483202.03*tUT1)+ 0.28*sind(228.2 + 960400.87*tUT1) - 0.28*sind(318.3 + 6003.18*tUT1) - 0.17*sind(217.6 - 407332.20*tUT1);
p = 0.9508 + 0.0518*cosd(134.9 + 477198.85*tUT1) + 0.0095*cosd(259.2 - 413335.38*tUT1) + 0.0078*cosd(235.7 + 890534.23*tUT1) + 0.0028*cosd(269.9 + 954397.70*tUT1);
e = 23.439291 - 0.0130032*tUT1 - (1.64e-7)*tUT1^2 + (5.04e-7)*tUT1^3 ;
rM = 6300.0/sind(p);
uVec = [cosd(phiEC)*cosd(lEC) cosd(e)*cosd(phiEC)*sind(lEC)-sind(e)*sind(phiEC) sind(e)*cosd(phiEC)*sind(lEC)+cosd(e)*sind(phiEC)];
rMOON = rM*(uVec);
