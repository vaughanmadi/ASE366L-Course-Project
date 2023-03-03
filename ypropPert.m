function [ydot] = ypropPert(t,y,V)
JD = 2458910;
muE = 398600.4415;
rE = 6378.1363;
wE = [0 0 2*pi/86164];

J2 = V(1);
J3 = V(2); % these load in my array of variables from V so that in ODE45, theres only one variable 'V' to propogate
Cr = V(4);
Cd = V(3);
Am = V(5);
muS = V(6);
muM = V(7);
t = JD + t./86400; % converting timespan into day fraction for later calculations
rE3 = sunpos(JD).*1.496e8; % CONVERT TO KM - dist from earth to sun add digits!!!
rEM = moonpos(JD); % dist from earth to moon
rES = [y(1) y(2) y(3)]; % dist from earth to spacecraft, given in problem
rS3 = rE3 - rES; % dist from s/c to sun, plug into thirdbod pert
rSM = rEM - rES; % dist from s/c to moon, plug into thirdbod pert

rmag = (y(1)^2 + y(2)^2 +y(3)^2)^(1/2);
ydot = zeros(size(y));
vSat = [y(4) y(5) y(6)];
aSUN = thirdbod(rE3,rS3,muS); % due to SUN third body
aMOON = thirdbod(rEM,rSM,muM); % due to MOON third body
aTBi = aSUN(1) + aMOON(1); 
aTBj = aSUN(2) + aMOON(2);
aTBk = aSUN(3) + aMOON(3);



aP = J2J3accel(rES,J2,J3,rE,muE); % due to J2 and J3
aJi = aP(1); 
aJj = aP(2); 
aJk = aP(3);

aSRP = srPressure(rES,rE3, Cr, Am); %due to solar radiation pressure
aSRPi = aSRP(1);
aSRPj = aSRP(2);
aSRPk = aSRP(3);

aD = dragPert(rES,vSat,Am,Cd,wE); %due to drag
aDi = aD(1);   
aDj = aD(2); 
aDk = aD(3);

ydot(1) = y(4); 
ydot(2) = y(5);
ydot(3) = y(6);
% add em up ;)
ydot(4) = -(muE/(rmag^3))*y(1) + aTBi + aJi + aSRPi + aDi;
ydot(5) = -(muE/(rmag^3))*y(2) + aTBj + aJj + aSRPj + aDj;
ydot(6) = -(muE/(rmag^3))*y(3) + aTBk + aJk + aSRPk + aDk; 
end