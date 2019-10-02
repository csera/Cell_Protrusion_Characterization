R = 3;   % Protein radius; nm
r = 0.5; % Small solute radius; nm
Posm = 2E5; % Osmotic P from the small solutes; Pa

Dmax = r*3;
D = linspace(0,Dmax);  % Dist b/w protein & flat wall

%Convert Posm to pN/nm^2 (1 Pa = 1E-6 pN/nm^2)
Posm = Posm * 1E-6;

F = -Posm*pi*(((D-r).^2)-(R+r)^2);

%Note on signs: since D has been defined to be +, F<0 indicates attraction
plot(D,F)
xlabel('Distance (nm)')
ylabel('Force (pN)')