% Version: 2.1      2019-08-09

% Numerical calculation of membrane protrusion bending stress & energy
% for an array of torus tube radii
% Breaking protrusion into 3 parts:
%       a) Base: half of the inner surface of a torus
%       b) Body: cylinder
%       c) Cap: hemisphere
%
% Exports the z slices and tension profile as .csv files

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
ETot = zeros(length(rt),3);             %Total bending E's for each seg

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
        
        %The "stress" (G per unit membrane area) along the base
        Ga(rChoice,a) = Kb*0.5*(((1/rt(rChoice))+(1/rl(a)))^2);
        ETot(rChoice,1) = ETot(rChoice,1) + Ga(rChoice,a);
    end
    
    %-------------------------------
    % Calc body bending E
    %-------------------------------
    %Dist along protrusion for body (abs meas)
    bodyEnd = ceil(n - (R*n/l));

    for a=baseEnd+1:bodyEnd
        Ga(rChoice,a) = Kb*0.5*((1/R)^2);
        ETot(rChoice,2) = ETot(rChoice,2) + Ga(rChoice,a);
    end

    %-------------------------------
    % Calc cap bending E
    %-------------------------------
    %Dist along protrusion for cap (abs meas)
    capEnd = n-1;
    
    for a=bodyEnd+1:capEnd
        Ga(rChoice,a) = Kb*2/(R^2);
        ETot(rChoice,3) = ETot(rChoice,3) + Ga(rChoice,a);
    end
end

%-------------------------------
% Graph the results
%-------------------------------

% Per 2019-08-07 work (see lab notebook), mechanical tension = bending
% energy per unit area

figure(1)
hold on
for rChoice=1:length(rt)
    plot(z,Ga(rChoice,:),'DisplayName',"r_t = "+num2str(rt(rChoice))+" nm")
end
legend
xlabel('Distance along protrusion (nm)')
ylabel('Tension (kT/nm^2)')
hold off


figure(2)
%Make labels for calculations with each rt
seriesLabels = strings(length(rt),1);
for a=1:length(rt)
    seriesLabels(a) = "rt = "+num2str(rt(a));
end

hold on

subplot(1,2,1);
bar(ETot,'stacked');
set(gca,'XTickLabel',seriesLabels)
legend('Base','Body','Cap')
ylabel('Total Segment Force (kT/nm)')

subplot(1,2,2);
bar(ETot);
set(gca,'XTickLabel',seriesLabels)      %Apply labels to each bar group
legend('Base','Body','Cap')
ylabel('Total Segment Force (kT/nm)')

hold off