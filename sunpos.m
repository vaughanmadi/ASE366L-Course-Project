function [rSUN] = sunpos(JDut1)
tUT1 = (JDut1-2451545)./36525;
lMs = 280.46 + 36000.771.*tUT1;
Msun = 357.52772333 + 35999.05034.*tUT1;
lECL = lMs + 1.914666471.*sind(Msun) + 0.019994643.*sind(2*Msun);

rsun = 1.000140612 - 0.016708617.*cosd(Msun) - 0.000139589.*cosd(2*Msun);
E = 23.439291 - 0.0130042.*tUT1;

rMOD = [rsun.*cosd(lECL) rsun.*cosd(E).*sind(lECL) rsun.*sind(E).*sind(lECL)];
   % answer is in AU and MOD
TT = tUT1/(365.25*100); 
x = deg2rad(2306.2181/3600).*TT + 0.30188.*TT.^2 + 0.017998.*TT.^3;
y = deg2rad(2004.3109/3600).*TT - 0.42665.*TT.^2 - 0.041833.*TT.^3;
z = deg2rad(2306.2181/3600).*TT + 1.09468.*TT.^2 + 0.018203.*TT.^3;

rSUN = R3(x)*R2(-y)*R3(z)*rMOD';
rSUN = rSUN';% now in GCRF and AU
end