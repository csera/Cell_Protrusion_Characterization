%Code from class sans text
%NOTE: For some reason, this produces a wave packet of huge amplitude
    %when done with dt >= ~1E-5

%model params
D = 10; %um^2/s (diffusion constant)
tTot = 1E-3; %s
L = 2; %length of cell in um
N0 = 100; %# molecules

%sim params
dt = 1E-6; %s
dx = 0.01; %um
k = D/(dx^2); %jumping rate const; let be uniform for all cells
M = L/dx;
N = zeros(ceil(tTot/dt),M); %# molecs at each point at each time
N(1,M/2) = N0;
for t=2:tTot/dt %time iterations, start at 2nd row (1st already defined)
    N(t,1) = N(t-1,1) + N(t-1,1+1)*k*dt - N(t-1,1)*k*dt; %1st space cell
    
    for x=2:M-1 %middle space
       N(t,x) = N(t-1,x) + N(t-1,x-1)*k*dt + N(t-1,x+1)*k*dt... 
            -2*N(t-1,x)*k*dt;
    end
    
    N(t,M) = N(t-1,M) + N(t-1,M-1)*k*dt - N(t-1,M)*k*dt; %last space cell
end

cell = 1:M;
plot(cell,N(100,:)) %plot ditr at t=100
xlabel('Box #')
ylabel('# molecules')

hold on %plot next thing on same graph
plot(cell,N(400,:),'-r') %plot ditr at t=400 in red
    %choose this t since you see avg displacement 2x at 4x the time in
    %diffusion
legend('t = 100','t = 400')