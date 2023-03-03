clear
format long 
%% Task 1 LEO
a = 6763; % KM
e = 0.001;
i = deg2rad(50); % rad
W = 0; 
w = 0; 
nu = 0;
muE = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
Cd = 2.0;
Cr = 1.5;
Am = 0.01;
muS = 1.327e11;
muM = 3903;
tp = greg2JD(2020,3,1,12,0,0); % March 1 2020 12:00:00 UTC - leap year????
[r1,v1] = orbtocart(a,e,i,W,w,nu,muE,0);
timespan = 0:600:432000; % 10-min increment to five days, in seconds
y = [r1(1) r1(2) r1(3) v1(1) v1(2) v1(3)]; % initial position

% propogate 2-body
odeoptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-12); 
[T, Y2b] =  ode45(@yprop2b, timespan, y, odeoptions,muE);
r2b = Y2b(:,1:3);

% propogate j2 + 2-body
V1 = variables(J2,0,0,0,0,0,0);
[T, YJ2] =  ode45(@ypropPert, timespan, y, odeoptions,V1);
rj2 = YJ2(:,1:3);
rDj2 = zeros(length(rj2),1);
for i = 1:length(rj2)
    rDj2(i) = norm(rj2(i,1:3)-r2b(i,1:3)); % taking magnitude of difference
end

% propogate j3 + 2-body
V2 = variables(0,J3,0,0,0,0,0);
[T, YJ3] =  ode45(@ypropPert, timespan, y, odeoptions,V2);
rj3 = YJ3(:,1:3);
rDj3 = zeros(length(rj3),1);
for i = 1:length(rj3)
    rDj3(i) = norm(rj3(i,1:3)-r2b(i,1:3)); 
end

% propogate drag + 2-body
V3 = variables(0,0,Cd,0,Am,0,0);
[T, YCd] =  ode45(@ypropPert, timespan, y, odeoptions,V3);
rCd = YCd(:,1:3);
rDCd = zeros(length(rCd),1);
for i = 1:length(rCd)
    rDCd(i) = norm(rCd(i,1:3)-r2b(i,1:3));
end

% propogate sun + 2-body
V4 = variables(0,0,0,0,0,muS,0);
[T, YSun] =  ode45(@ypropPert, timespan, y, odeoptions,V4);
rSun = YSun(:,1:3);
rDSun = zeros(length(rSun),1);
for i = 1:length(rSun)
    rDSun(i) = norm(rSun(i,1:3)-r2b(i,1:3));
end

% propogate moon + 2-body
V5 = variables(0,0,0,0,0,0,muM);
[T, YMoon] =  ode45(@ypropPert, timespan, y, odeoptions,V5);
rMoon = YMoon(:,1:3);
rDMoon = zeros(length(rMoon),1);
for i = 1:length(rMoon)
    rDMoon(i) = norm(rMoon(i,1:3)-r2b(i,1:3));
end

% propogate SRP + 2-body
V6 = variables(0,0,0,Cr,Am,0,0);
[T, YSRP] =  ode45(@ypropPert, timespan, y, odeoptions,V6);
rSRP = YSRP(:,1:3);
rDSRP = zeros(length(rSRP),1);
for i = 1:length(rSRP)
    rDSRP(i) = norm(rSRP(i,1:3)-r2b(i,1:3));
end

days = 0:5/(length(T)-1):5; % create axis for graphing

figure(1)
set(0,'DefaultAxesFontName','Arial') %wtf font to use
set(0,'DefaultAxesFontSize',12)

semilogy(days,rDj2)
hold on
semilogy(days,rDj3)
semilogy(days,rDCd)
semilogy(days,rDSun)
semilogy(days,rDMoon)
semilogy(days,rDSRP)
hold off
legend ('J2','J3','Atmospheric Drag','Sun Third Body','Moon Third Body','Solar Radiation Pressure','location','southeast')
xlabel('Time from Epoch (days)')
ylabel('Difference (km)')
title('LEO')
%% Task 1 MEO
a = 26560; 
e = 0.001;
i = deg2rad(55); %rad
W = 0;
w = 0;

[r1,v1] = orbtocart(a,e,i,W,w,nu,muE,0);
timespan = 0:600:432000; % 10-min increment to five days, in seconds
y = [r1(1) r1(2) r1(3) v1(1) v1(2) v1(3)]; % initial position

% propogate 2-body
odeoptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-12); 
[T, Y2b] =  ode45(@yprop2b, timespan, y, odeoptions,muE);
r2b = Y2b(:,1:3);

% propogate j2 + 2-body
V1 = variables(J2,0,0,0,0,0,0);
[T, YJ2] =  ode45(@ypropPert, timespan, y, odeoptions,V1);
rj2 = YJ2(:,1:3);
rDj2 = zeros(length(rj2),1);
for i = 1:length(rj2)
    rDj2(i) = norm(rj2(i,1:3)-r2b(i,1:3)); % taking magnitude of difference
end

% propogate j3 + 2-body
V2 = variables(0,J3,0,0,0,0,0);
[T, YJ3] =  ode45(@ypropPert, timespan, y, odeoptions,V2);
rj3 = YJ3(:,1:3);
rDj3 = zeros(length(rj3),1);
for i = 1:length(rj3)
    rDj3(i) = norm(rj3(i,1:3)-r2b(i,1:3)); 
end

% propogate drag + 2-body
V3 = variables(0,0,Cd,0,Am,0,0);
[T, YCd] =  ode45(@ypropPert, timespan, y, odeoptions,V3);
rCd = YCd(:,1:3);
rDCd = zeros(length(rCd),1);
for i = 1:length(rCd)
    rDCd(i) = norm(rCd(i,1:3)-r2b(i,1:3));
end

% propogate sun + 2-body
V4 = variables(0,0,0,0,0,muS,0);
[T, YSun] =  ode45(@ypropPert, timespan, y, odeoptions,V4);
rSun = YSun(:,1:3);
rDSun = zeros(length(rSun),1);
for i = 1:length(rSun)
    rDSun(i) = norm(rSun(i,1:3)-r2b(i,1:3));
end

% propogate moon + 2-body
V5 = variables(0,0,0,0,0,0,muM);
[T, YMoon] =  ode45(@ypropPert, timespan, y, odeoptions,V5);
rMoon = YMoon(:,1:3);
rDMoon = zeros(length(rMoon),1);
for i = 1:length(rMoon)
    rDMoon(i) = norm(rMoon(i,1:3)-r2b(i,1:3));
end

% propogate SRP + 2-body
V6 = variables(0,0,0,Cr,Am,0,0);
[T, YSRP] =  ode45(@ypropPert, timespan, y, odeoptions,V6);
rSRP = YSRP(:,1:3);
rDSRP = zeros(length(rSRP),1);
for i = 1:length(rSRP)
    rDSRP(i) = norm(rSRP(i,1:3)-r2b(i,1:3));
end

days = 0:5/(length(T)-1):5; % create axis for graphing

figure(2)
set(0,'DefaultAxesFontName','Arial') 
set(0,'DefaultAxesFontSize',12)

semilogy(days,rDj2)
hold on
semilogy(days,rDj3)
semilogy(days,rDCd)
semilogy(days,rDSun)
semilogy(days,rDMoon)
semilogy(days,rDSRP)
hold off
legend ('J2','J3','Atmospheric Drag','Sun Third Body','Moon Third Body','Solar Radiation Pressure','location','southeast')
xlabel('Time from Epoch (days)')
ylabel('Difference (km)')
title('MEO')
%% Task 1 GEO
a = 42164;
e = 0.01;
i = deg2rad(0.5); %rad
W = deg2rad(-120); % rad
w = 0;
[r1,v1] = orbtocart(a,e,i,W,w,nu,muE,0);
timespan = 0:600:432000; % 10-min increment to five days, in seconds
y = [r1(1) r1(2) r1(3) v1(1) v1(2) v1(3)]; % initial position

% propogate 2-body
odeoptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-12); 
[T, Y2b] =  ode45(@yprop2b, timespan, y, odeoptions,muE);
r2b = Y2b(:,1:3);

% propogate j2 + 2-body
V1 = variables(J2,0,0,0,0,0,0);
[T, YJ2] =  ode45(@ypropPert, timespan, y, odeoptions,V1);
rj2 = YJ2(:,1:3);
rDj2 = zeros(length(rj2),1);
for i = 1:length(rj2)
    rDj2(i) = norm(rj2(i,1:3)-r2b(i,1:3)); % taking magnitude of difference
end

% propogate j3 + 2-body
V2 = variables(0,J3,0,0,0,0,0);
[T, YJ3] =  ode45(@ypropPert, timespan, y, odeoptions,V2);
rj3 = YJ3(:,1:3);
rDj3 = zeros(length(rj3),1);
for i = 1:length(rj3)
    rDj3(i) = norm(rj3(i,1:3)-r2b(i,1:3)); 
end

% propogate drag + 2-body
V3 = variables(0,0,Cd,0,Am,0,0);
[T, YCd] =  ode45(@ypropPert, timespan, y, odeoptions,V3);
rCd = YCd(:,1:3);
rDCd = zeros(length(rCd),1);
for i = 1:length(rCd)
    rDCd(i) = norm(rCd(i,1:3)-r2b(i,1:3));
end

% propogate sun + 2-body
V4 = variables(0,0,0,0,0,muS,0);
[T, YSun] =  ode45(@ypropPert, timespan, y, odeoptions,V4);
rSun = YSun(:,1:3);
rDSun = zeros(length(rSun),1);
for i = 1:length(rSun)
    rDSun(i) = norm(rSun(i,1:3)-r2b(i,1:3));
end

% propogate moon + 2-body
V5 = variables(0,0,0,0,0,0,muM);
[T, YMoon] =  ode45(@ypropPert, timespan, y, odeoptions,V5);
rMoon = YMoon(:,1:3);
rDMoon = zeros(length(rMoon),1);
for i = 1:length(rMoon)
    rDMoon(i) = norm(rMoon(i,1:3)-r2b(i,1:3));
end

% propogate SRP + 2-body
V6 = variables(0,0,0,Cr,Am,0,0);
[T, YSRP] =  ode45(@ypropPert, timespan, y, odeoptions,V6);
rSRP = YSRP(:,1:3);
rDSRP = zeros(length(rSRP),1);
for i = 1:length(rSRP)
    rDSRP(i) = norm(rSRP(i,1:3)-r2b(i,1:3));
end

days = 0:5/(length(T)-1):5; % create axis for graphing

figure(3)
set(0,'DefaultAxesFontName','Arial') %wtf font to use
set(0,'DefaultAxesFontSize',12)

semilogy(days,rDj2)
hold on
semilogy(days,rDj3)
semilogy(days,rDCd)
semilogy(days,rDSun)
semilogy(days,rDMoon)
semilogy(days,rDSRP)
hold off
legend ('J2','J3','Atmospheric Drag','Sun Third Body','Moon Third Body','Solar Radiation Pressure','location','southeast')
xlabel('Time from Epoch (days)')
ylabel('Difference (km)')
title('GEO')
%% Task 1 Molniya
a = 26000;
e = 0.72;
i = deg2rad(75); %rad
W = deg2rad(90); %rad
w = deg2rad(-90); %rad
[r1,v1] = orbtocart(a,e,i,W,w,nu,muE,0);
timespan = 0:600:432000; % 10-min increment to five days, in seconds
y = [r1(1) r1(2) r1(3) v1(1) v1(2) v1(3)]; % initial position

% propogate 2-body
odeoptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-12); 
[T, Y2b] =  ode45(@yprop2b, timespan, y, odeoptions,muE);
r2b = Y2b(:,1:3);

% propogate j2 + 2-body
V1 = variables(J2,0,0,0,0,0,0);
[T, YJ2] =  ode45(@ypropPert, timespan, y, odeoptions,V1);
rj2 = YJ2(:,1:3);
rDj2 = zeros(length(rj2),1);
for i = 1:length(rj2)
    rDj2(i) = norm(rj2(i,1:3)-r2b(i,1:3)); % taking magnitude of difference
end

% propogate j3 + 2-body
V2 = variables(0,J3,0,0,0,0,0);
[T, YJ3] =  ode45(@ypropPert, timespan, y, odeoptions,V2);
rj3 = YJ3(:,1:3);
rDj3 = zeros(length(rj3),1);
for i = 1:length(rj3)
    rDj3(i) = norm(rj3(i,1:3)-r2b(i,1:3)); 
end

% propogate drag + 2-body
V3 = variables(0,0,Cd,0,Am,0,0);
[T, YCd] =  ode45(@ypropPert, timespan, y, odeoptions,V3);
rCd = YCd(:,1:3);
rDCd = zeros(length(rCd),1);
for i = 1:length(rCd)
    rDCd(i) = norm(rCd(i,1:3)-r2b(i,1:3));
end

% propogate sun + 2-body
V4 = variables(0,0,0,0,0,muS,0);
[T, YSun] =  ode45(@ypropPert, timespan, y, odeoptions,V4);
rSun = YSun(:,1:3);
rDSun = zeros(length(rSun),1);
for i = 1:length(rSun)
    rDSun(i) = norm(rSun(i,1:3)-r2b(i,1:3));
end

% propogate moon + 2-body
V5 = variables(0,0,0,0,0,0,muM);
[T, YMoon] =  ode45(@ypropPert, timespan, y, odeoptions,V5);
rMoon = YMoon(:,1:3);
rDMoon = zeros(length(rMoon),1);
for i = 1:length(rMoon)
    rDMoon(i) = norm(rMoon(i,1:3)-r2b(i,1:3));
end

% propogate SRP + 2-body
V6 = variables(0,0,0,Cr,Am,0,0);
[T, YSRP] =  ode45(@ypropPert, timespan, y, odeoptions,V6);
rSRP = YSRP(:,1:3);
rDSRP = zeros(length(rSRP),1);
for i = 1:length(rSRP)
    rDSRP(i) = norm(rSRP(i,1:3)-r2b(i,1:3));
end

days = 0:5/(length(T)-1):5; % create axis for graphing

figure(4)
set(0,'DefaultAxesFontName','Arial') 
set(0,'DefaultAxesFontSize',12)

semilogy(days,rDj2)
hold on
semilogy(days,rDj3)
semilogy(days,rDCd)
semilogy(days,rDSun)
semilogy(days,rDMoon)
semilogy(days,rDSRP)
hold off
legend ('J2','J3','Atmospheric Drag','Sun Third Body','Moon Third Body','Solar Radiation Pressure','location','southeast')
xlabel('Time from Epoch (days)')
ylabel('Difference (km)')
title('Molniya')
%% Task 2 
J2 = 0.0010826267;
J3 = -0.0000025327;
Cd = 2.0;
Cr = 1.5;
Am = 0.01;
muS = 1.327e11;
muM = 3903;
JD = 2458910; % same start time task 1
lat1 = deg2rad(30.30); %N, deg
long1 = deg2rad(-120.6); %W, deg
lat2 = deg2rad(30.3); %N, deg
long2 = deg2rad(-97.7); %W, deg
norbits = 10;
nu1 = deg2rad(20.00); %deg, angle between perigee and pos. burnout
a = 6500; %km
e = 0.001;
% theta = lat, lambda = long
wE = 2*pi/86164; % rad/s
muE = 398600.4415; 
rE = 6378.1363; % km
g = 9.807*0.001; %m/s^2 -> km/s^2
n = sqrt(muE/(a^3));
p = a*(1-e^2);

% initial values
T = 2*pi*(1/n);
E = 2*atan2(sqrt(1-e)*tan(nu1/2),sqrt(1+e));
t1 = (T/(2*pi))*(E - e*sin(E)); % time since periapsis at burnout
long2e = long2 + norbits*wE*T;


% 'first iteration'
t2et1 = (T/(2*pi))*(long2e-long1); % eq 21 initial guess for change in time
deltalong2e1 = long2e - long1 + wE*t2et1;
nu2e = nu1 + acos(sin(lat2)*sin(lat1)+cos(lat2)*cos(lat1)*cos(deltalong2e1));
Enew = 2*atan2(sqrt(1-e)*tan(nu2e/2),sqrt(1+e));
t2e = (T/(2*pi))*(Enew - e*sin(Enew)); % eq 8 new time2
updatedt2et1 = t2e - t1;

az = asin((sin(deltalong2e1)*cos(lat2))/sin(nu2e-nu1));
i = acos(sin(az)*cos(lat1));

C1 = -((3/2)*J2*sqrt(muE))/(rE^(3/2));
dRAAN = C1*((rE/p)^2)*((rE/a)^(3/2))*cos(i);
dt = norbits*T + updatedt2et1; 
deltaRAAN = dRAAN*dt; % oblateness correction

dw = (3*n*(rE^2)*J2)*(4-5*(sin(i)^2))/(4*p^2);
deltaw = dw*dt; % oblateness correction
w = asin(sin(lat1)/sin(i)) - nu1;
deltalat = (sin(i)*cos(w+nu2e)*deltaw)/cos(lat2); % corrected
deltalong = (cos(i)*(sec(w+nu2e)^2)*(deltaw))/(1+(cos(i)^2)*(tan(w+nu2e)^2)) + deltaRAAN; % updated
lat2correct = lat2 - deltalat;
long2correct = long2 - deltalong;


long2e = long2correct + norbits*wE*T;
for j = 1:3 %for 4 iterations
deltalong2e1 = long2e - long1 + wE*updatedt2et1;
nu2e = nu1 + acos(sin(lat2)*sin(lat1)+cos(lat2)*cos(lat1)*cos(deltalong2e1));
Enew = 2*atan2(sqrt(1-e)*tan(nu2e/2),sqrt(1+e));
t2e = (T/(2*pi))*(Enew - e*sin(Enew)); 
az = asin((sin(deltalong2e1)*cos(lat2))/sin(nu2e-nu1));
i = acos(sin(az)*cos(lat1));
updatedt2et1 = t2e - t1; 
w = asin(sin(lat1)/sin(i)) - nu1;
end 

%final calcs for finding sat pos and RAAN
deltalongN1 = atan2(sin(lat1)*sin(az),cos(az));
longNref = long1 - deltalongN1;
tGMST = 0; % assumed zero at burnout
RAAN = tGMST + longNref;

[r1,v1] = orbtocart(a,e,i,RAAN+90,w,nu2e,muE,0);
% v1 is the answer for 2.1

RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);
y1 = [r1(1) r1(2) r1(3) v1(1) v1(2) v1(3)];

timespan = 0:30:norbits*T;
odeoptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);


[time,yprop] = ode45(@ypropJ2,timespan,y1,odeoptions,J2);
rITRF = yprop(:,1:3);
rGCRF = zeros(length(rITRF),3);


jdVec = JD + timespan./86400;
geolat = zeros(length(rITRF),1);
geolong = zeros(length(rITRF),1);

for j = 1:length(rITRF)
thetaG = wE*(time(j));
rGCRF(j,1:3) = R3(thetaG)*rITRF(j,1:3)'; % convert ITRF to GCRF
geolat(j) = rad2deg(asin(rGCRF(j,3)/norm(rGCRF(j,1:3))));
geolong(j) = rad2deg(atan2(rGCRF(j,2),rGCRF(j,1)));
end
long1 = rad2deg(long1); % adjusting longitude to fit groundtrack map
lat1 = rad2deg(lat1);
long2 = rad2deg(long2);
lat2 = rad2deg(lat2);

% plot the groundtrack
load earth_coastline.mat 
figure(5)
set(0,'DefaultAxesFontName','Arial') %wtf font to use
set(0,'DefaultAxesFontSize',12)
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
hold on
plot(geolong,geolat,'b.','MarkerSize',8)
plot(long1,lat1,'r.','MarkerSize',20) 
plot(long2,lat2,'g.','MarkerSize',20)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title('Ground Track of Two Body Orbit with J2')
hold off
legend('Coastline','Groundtrack','Burnout Position','Observation Point','location','Southeast')
axis equal
xlim([-180,180])
ylim([-90,90])

% use high fidelity propagator
V = variables(J2,J3,Cd,Cr,Am,muS,muM);
[time,ypropHF] = ode45(@ypropPert,timespan,y1,odeoptions,V);
rITRFHF = ypropHF(:,1:3);
rGCRFHF = zeros(length(rITRFHF),3);

geolatHF = zeros(length(rITRFHF),1);
geolongHF = zeros(length(rITRFHF),1);

for j = 1:length(rITRFHF)
thetaG = wE*(time);
rGCRFHF(j,1:3) = R3(thetaG(j))*rITRFHF(j,1:3)'; % convert ITRF to GCRF
geolatHF(j) = rad2deg(asin(rGCRFHF(j,3)/norm(rGCRFHF(j,1:3))));
geolongHF(j) = rad2deg(atan2(rGCRFHF(j,2),rGCRFHF(j,1)));
end

figure(6)
set(0,'DefaultAxesFontName','Arial') %wtf font to use
set(0,'DefaultAxesFontSize',12)
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
hold on
plot(geolong,geolat,'b.','MarkerSize',8)
plot(long1,lat1,'r.','MarkerSize',20) 
plot(long2,lat2,'g.','MarkerSize',20)
plot(geolongHF,geolatHF,'c.','MarkerSize',8)
hold off
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title('Ground Track of Two Body Orbit with all Special Pertubations')
legend('Coastline','Groundtrack','Burnout Position','Observation Point','High Fidelity Propagation','location','Southeast')
axis equal
xlim([-180,180])
ylim([-90,90])
% why does it look like this </3