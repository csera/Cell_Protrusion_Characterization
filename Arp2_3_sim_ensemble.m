% SIMULATING ARP2/3 ABM TRAJECTORIES
% Curtis Sera, Welch Lab
% v1.0, 2020-01-16
%
% This code is adapted from Arp2_3_sim.m v1.1
% 
% Motivation: Bp utilizes the Ena/VASP machinery for ABM and thus produces 
% straight tails.  However, Bt utilizes the branching Arp2/3 machinery, and
% this has been empirically found to yield curved tails.  Why is this?
%   This seemed intuitive to me: branch --> curve. Matt didn't find this so
%   intuitive though, and upon closer inspection, it certainly could seem
%   surprising that the Arp2/3 tails are neatly curved rather than, say,
%   zig-zagged.
%
% Program overview:
% Like Arp2_3_sim.m, this script will sim ABM construction as a Poisson
% process in which null events are addition of actin monomers to the tail
% (ie tail lengthens in whatever dxn it's currently pointing) and point
% events are binding of Arp2/3 (ie tail growht direction deviates 70
% degrees).  Whether the tail deviates "left" or "right" will be a 50/50
% coin toss.
% This script will seek to examine ensemble behavior rather than that of
% individual filaments, and it will thus use larger nSims to reach more
% representative RMS and mean paths.
%   These will be plotted against each other for each PArp without the
%   clutter of all the individual paths

clear
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
nSims = 250;    %# of runs to simulate
    % Bt width/f-actin diameter ~= 500 nm/2 nm = 250 f-actin

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

figure(1)
hold on
for p=1:pRuns
    plot(rmsX(:,p),rmsY(:,p),'k','Color',colors{p},'LineWidth',2,...
        'DisplayName',"RMS path, P = "+num2str(PArp(p)));
end
title("RMS paths")
lgd = legend();
lgd.Location = 'northwest';
hold off

figure(2)
hold on
for p=1:pRuns
    plot(mX(:,p),mY(:,p),':k','Color',colors{p},'LineWidth',2,...
        'DisplayName',"Mean path, P = "+num2str(PArp(p)));
end
title("Mean paths")
lgd = legend();
lgd.Location = 'northwest';
hold off