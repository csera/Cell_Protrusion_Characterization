% SIMULATING ARP2/3 ABM TRAJECTORIES
% Curtis Sera, Welch Lab
% v1.1, 2020-01-15

% Motivation: Bp utilizes the Ena/VASP machinery for ABM and thus produces 
% straight tails.  However, Bt utilizes the branching Arp2/3 machinery, and
% this has been empirically found to yield curved tails.  Why is this?
%   This seemed intuitive to me: branch --> curve. Matt didn't find this so
%   intuitive though, and upon closer inspection, it certainly could seem
%   surprising that the Arp2/3 tails are neatly curved rather than, say,
%   zig-zagged.
%
% Program overview:
% This script will simulate ABM construction as a Poisson process in which
% null events are addtions of actin monomers to the tail (ie tail
% lengthens) and point events are binding of Arp2/3 (ie the tail redirects
% its growth to be 70 degrees deviated).  There will be chance 0.5 that
% deviation is to the "left" and chance 0.5 that it is to the "right".
% This first-pass will only look at nSims runs of an isolated filament in a
% 2D plane

clear all
close all

% Assumed model params
PArp = [0.001; 0.01; 0.1; 0.3;0.8;0.95];
    %P that Arp2/3 is added in rather than a new actin monomer
pRuns = size(PArp,1);
PRt = 0.5;      %P that Arp2/3 addition makes turn "right"
LAct = 1;       %length added by addition of actin monomer
ArpAng = 70 * pi/180;       %Angle deviation in rad caused by Arp2/3

% Sim settings
tEnd = 1000;    %# time steps
nSims = 20;      %# of runs to simulate

% Sim vars
a = zeros(tEnd,nSims,pRuns);      %Angle of trajectory
x = zeros(tEnd,nSims,pRuns);      %x coordinate
y = zeros(tEnd,nSims,pRuns);      %y coordinate
arp = rand(tEnd,nSims,pRuns);     %Arp2/3 added if arp(t} <= PArp

rmsX = zeros(tEnd,pRuns);
rmsY = zeros(tEnd,pRuns);
mX = zeros(tEnd,pRuns);
mY = zeros(tEnd,pRuns);

%For pRun, iterate at each time point through all the sims to get new
%positions
for p=1:pRuns
    for t=2:tEnd 
        for n=1:nSims
            %First adjust elongation angle as appropriate
            if arp(t,n,p)<= PArp(p)
                RL = rand;
                if RL<0.5
                    a(t,n,p) = a(t-1,n,p)+ArpAng;
                else
                    a(t,n,p) = a(t-1,n,p)-ArpAng;
                end
            else
                a(t,n,p) = a(t-1,n,p);
            end

            x(t,n,p) = x(t-1,n,p) + LAct*cos(a(t,n,p) * pi/180);
            y(t,n,p) = y(t-1,n,p) + LAct*sin(a(t,n,p) * pi/180);

        end
        
        %Get RMS of this time point
        rmsX(t,p) = sqrt(mean(x(t,:,p).^2));
        rmsY(t,p) = sqrt(mean(y(t,:,p).^2));
        
        %Get mean of this time point
        mX(t,p) = mean(x(t,:,p));
        mY(t,p) = mean(y(t,:,p));
        
    end
end

%Plot results
colors = {'r','m','y','g','c','b'};

for p=1:pRuns
    figure(p)
    hold on
    simsPlot = plot(x(:,:,p),y(:,:,p),'Color',colors{p});
    rmsSimsPlot = plot(rmsX(:,p),rmsY(:,p),'k','LineWidth',2,...
        'DisplayName','RMS position');
    mSimsPlot = plot(mX(:,p),mY(:,p),':k','LineWidth',2,...
        'DisplayName','Mean position');
    title("PArp = "+PArp(p))
    lgd = legend([rmsSimsPlot,mSimsPlot]);
    lgd.Location = 'southwest';
end

