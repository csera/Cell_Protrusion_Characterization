% Taking "discrete" derivatives
% Curtis Sera, Welch Lab
% 2020-01-17, v1.0
%
% I made some code for Arp2_3_sim_ensemble.m v1.2 that attempts to take the
% 1st and 2nd derivative of angle data over time.  However, the angle data
% that I analyzed with this code was pretty comlex because I generated it
% using a 2 layer stochastic process for a long series.  That is, I have no
% easy way of verifying that my code is doing what it's supposed to.
% This script will check my methodology by generating data from simple
% functions which have easily identifiable derivatives.

close all

% -----------------------------------------
% Generate the data
% -----------------------------------------
step = 0.1;
x = -5:step:5;
stop = size(x,2);

C = x.^3;
Ex = exp(2*x);
S = sin(x);

% -----------------------------------------
% Apply my code for taking discrete derivatives
% -----------------------------------------
dC = discDeriv(C,step,1);
ddC = discDeriv(dC,step,2);

dEx = discDeriv(Ex,step,1);
ddEx = discDeriv(dEx,step,2);

dS = discDeriv(S,step,1);
ddS = discDeriv(dS,step,2);

% -----------------------------------------
% Plot
% -----------------------------------------
% Graphing the origianl fxn, my computed derivatives, and comparing those
% to the actual derivative functions.
% Note: Plots of odd derivatives will all have their curves shifted left by
% 1 step due to how I managed the indexing.  This is expected.
figure(1)
hold on
plot(x,C,'DisplayName','C (x^3)')
plot(x,dC,'DisplayName','dC')
plot(x,3*x.^2,'DisplayName','3x^2')
plot(x,ddC,'DisplayName','ddC')
plot(x,6*x,'DisplayName','6x')
legend

figure(2)
hold on
plot(x,Ex,'DisplayName','Ex')
plot(x,dEx,'DisplayName','dEx')
plot(x,2*exp(2*x),'DisplayName','2e^(2x)')
plot(x,ddEx,'DisplayName','ddEx')
plot(x,4*exp(2*x),'DisplayName','4e^(2x)')
legend

figure(3)
hold on
plot(x,S,'DisplayName','C')
plot(x,dS,'DisplayName','dC')
plot(x,cos(x),'DisplayName','cos(x)')
plot(x,ddS,'DisplayName','ddC')
plot(x,-sin(x),'DisplayName','-sin(x)')
legend

function d = discDeriv(a,step,odd)
    stop = size(a,2);
    d = zeros(1,stop);
    
    if odd==1
        fprintf('odd \n')
        %Find odd derivative. Shift times back by half a step
        for n=1:stop-1
            d(n) = (a(n+1) - a(n))/step;
        end
        
    else
        fprintf('else \n')
        %Find even derivative. Times from dmA already match this
        for n=2:stop-1
            d(n) = (a(n) - a(n-1))/step;
        end
    end

end
