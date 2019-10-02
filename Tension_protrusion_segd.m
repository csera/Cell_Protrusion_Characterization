% Version: 1.0      2019-08-07

% Based on Gbend_protrusion2.m v2.0.1
% Numerical calculation of protrusion membrane tension profile based on
% the bending energy profile according to Lipowsky's (2014) relation.
% Protrusion segmented in 3 parts:
%       a) Base: half of the inner surface of a torus
%       b) Body: cylinder
%       c) Cap: hemisphere

clear;

%-------------------------------
%   Define system props
%-------------------------------
%Define membrane properties
Kb = 15;                %Bending modulus; kT
thickness = 5;          %Membrane thickness; nm

%Define projection properties
%Universal props
R = 250;                %Radius of the projection; nm
l = 5000;               %Protrusion length; nm
%Base/torus props
rt = [50,100,200,400,600];  %Radius of the torus tube; nm
%Body props
lCyl = l-(2*rt);        %Protrusion body length (ie w/o cap or base); nm
%Note: cap is hemisphere of radius R => no unique props

%-------------------------------
% Set up key vars
%-------------------------------
%Make our line for moving along the protrusion length
n = l*2;                %number of bins to split line into
z = linspace(0,l,n);    %array of bins; bin step size l/n nm; units: nm
    z = z(2:n);         %Remove 0 bin to avoid inf problems
    %z will be right-side inclusive. Ie bin 2 will cover (l/n,2*l/n]
Ga = zeros(length(rt),length(z));       %Bending energy/area; kT/nm^2
drl = zeros(length(rt),length(z));      %d(rl)/dz

%-------------------------------
% Calc bending E's
%-------------------------------
%Iterate through ri and calculate the stress profile using each ri val
for rChoice=1:length(rt)
    %-------------------------------
    % Calc base bending E
    %-------------------------------
    %Index of z where base ends for ri(rChoice)
    baseEnd = ceil(rt(rChoice)*n/l);
    
    rl = zeros(1,baseEnd);      %Dist b/w protrusion center & sides; nm
        %Note that this is a fxn of z
    
    for a=1:baseEnd
        rl(a) = R + rt(rChoice)...
            - (rt(rChoice)*cos(asin((rt(rChoice)-z(a))/rt(rChoice))));
        drl(rChoice,a) = -((rt(rChoice)-z(a))/rt(rChoice))...
            *((1-((rt(rChoice)-z(a))/rt(rChoice))^2)^(-0.5));
        
        %The "stress" (G per unit membrane area) along the base
        Ga(rChoice,a) = Kb*0.5*(((1/rt(rChoice))+(1/rl(a)))^2);
    end
    
    %-------------------------------
    % Calc body bending E
    %-------------------------------
    %Dist along protrusion for body (abs meas)
    bodyEnd = ceil(n - (R*n/l));

    for a=baseEnd+1:bodyEnd
        Ga(rChoice,a) = Kb*0.5*((1/R)^2);
    end

    %-------------------------------
    % Calc cap bending E
    %-------------------------------
    %Dist along protrusion for cap (abs meas)
    capEnd = n-1;
    
    for a=bodyEnd+1:capEnd
        Ga(rChoice,a) = Kb*2/(R^2);
    end
end

%-------------------------------
% Transform results 
%-------------------------------


%-------------------------------
% Graph the results
%-------------------------------
stress = Ga./thickness;  %Stress (E/V) per strip; kT/(nm^3)

figure(1)
hold on
baseEnd = ceil(max(rt)*n/l)*1.5;
for rChoice=1:length(rt)
    plot(z(1:baseEnd),drl(rChoice,1:baseEnd),'DisplayName',"r_t = "+...
        num2str(rt(rChoice))+" nm")
end
legend
xlabel('Distance along protrusion (nm)')
ylabel('d')
hold off

hold off
