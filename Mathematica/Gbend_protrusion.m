% Version: 1.0      2019-07-09

% Numerical calculation of membrane protrusion bending stress & energy
% for a single, chosen torus tube radius
% Breaking protrusion into 3 parts:
%       a) Base: half of the inner surface of a torus
%       b) Body: cylinder
%       c) Cap: hemisphere

clear;

%-------------------------------
%   Define system props
%-------------------------------
%Define membrane properties
Kb = 15;            %Bending modulus; kT
thickness = 5;      %Membrane thickness; nm

%Define projection properties
%Universal props
R = 250;            %Radius of the projection; nm
l = 5000;           %Protrusion length; nm
%Base/torus props
ri = 200;           %Radius of the torus tube; nm
rt = R+ri;        %Dist from torus center to tube center; nm
%Body props
lCyl = l-(2*ri);    %Protrusion body length (ie w/o cap or base); nm
%Note: cap is hemisphere of radius R => no unique props

%Make our line for moving along the protrusion length
n = l*2;                %number of bins to split line into
dz = l/n;               %dist bins step over; nm
z = linspace(0,l,n);    %array of bins; bin step size l/n nm; units: nm
    z = z(2:n);         %Remove 0 bin to avoid inf problems
    %z will be right-side inclusive. Ie bin 2 will cover (l/n,2*l/n]

%-------------------------------
% Calc base bending E
%-------------------------------
%Dist along protrusion for base (absolute meas)
zBase = z(1:ceil(ri*n/l));    %units: nm

%Dist b/w protrusion center & sides as fxn of z:
rl = R + ri - (ri*cos(asin((ri-zBase)/ri)));
    %Check: yes, rl vs z plots a circle
drl = -((ri-zBase)/ri).*((1-((ri-zBase)/ri).^2).^(-0.5)); %d(rl)/dz
    %exp checks out

%The "stress" (G per unit membrane area) along the base
sBase = Kb*0.5*(((1/ri)+(1./rl)).^2);
sBaseStrip = Kb*0.5*(((1/ri)+(1./rl)).^2)...
    *pi.*rl.*sqrt(1+drl.^2)*dz;

% figure(1)
% plot(zBase,sBase,'-r')
% xlabel("z (nm)")
% ylabel("G per unit area (kT)")
% 
% figure(2)
% plot(zBase,sBaseStrip,'-b')
% xlabel("z (nm)")
% ylabel("G per strip (kT)")

%The totaled bending energy of the base
sBaseSum = sum(sBaseStrip);

%-------------------------------
% Calc body bending E
%-------------------------------
%Dist along protrusion for body (abs meas)
zBody = z(ceil(ri*n/l)+1:ceil(n - (R*n/l)));

%The "stress" along the body (note: cylinder => one r is inf)
sBody = zeros(1,length(zBody));
sBodyStrip = sBody;

for a=1:length(zBody)
    sBody(a) = Kb*0.5*((1/R)^2);
    sBodyStrip(a) = Kb*pi*dz/R;
end

sBodySum = sum(sBody);
GBodyTheor = Kb*pi*(l-2*ri)/R;  %Theoretical eq for G of whole body

%-------------------------------
% Calc cap bending E
%-------------------------------
%Dist along protrusion for cap (abs meas)
zCap = z(ceil(n - (R*n/l))+1:n-1);

sCap = zeros(1,length(zCap));
%sCapStrip = sCap;              %I'm going to ignore this for now

for b=1:length(zCap)
    sCap(b) = Kb*2/(R^2);
end

GCapTheor = 4*pi*Kb;

%-------------------------------
% Graph the results
%-------------------------------
sFull = [sBase,sBody,sCap];     %G_bend per unit area per strip; kT/(nm^2)
stressFull = sFull./thickness;  %Stress (E/V) per strip; kT/(nm^3)

figure(1)
plot(z,sFull,'-b')
xlabel('z (nm)')
ylabel('G per unit area (kT/nm^2)')

figure(2)
plot(z,stressFull,'-r')
xlabel('z (nm)')
ylabel('Stress (kT/nm^3)')
