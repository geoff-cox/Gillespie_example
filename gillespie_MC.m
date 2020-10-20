function gillespie_MC()

close all;

N_sims = 100;   % number of simulations

t0 = 0;                     % initial time
tF = 30;                    % final time
dt_avg = 0.0455;            % average gillespie time step (computed separately)
time_std = t0:dt_avg:tF;    % time partition
N_t = length(time_std);     % number of time steps

gamma = 0.1;   % recovery rate
beta = 1.3;    % infection rate

N = 350;    % initial population
I0 = 6;     % initial infected number

% solve system of odes for given parameters using ode45
[t_ode,soln] = ode45(@(t,y) odefcn(y,gamma,beta), [t0,tF], [N-I0 I0 0]);
S_ode = soln(:,1);
I_ode = soln(:,2);
R_ode = soln(:,3);

% initialize stochastic group evolution vectors
S_avg = zeros(N_t,1);        % susceptible
I_avg = zeros(N_t,1);        % infected
R_avg = zeros(N_t,1);        % recovered
c  = zeros(N_t,1);           % frequency counter

% initial values
S_avg(1) = N-I0;
I_avg(1) = I0;
R_avg(1) = 0;

% -------------------- simulation loop --------------------
for sim = 1:N_sims
    t = 0;
    % initialize stochastic group evolution vectors
    S = zeros(N_t,1);        % susceptible
    I = zeros(N_t,1);        % infected
    R = zeros(N_t,1);        % recovered
    time = zeros(N_t,1);     % time
    
    % initialize time = 0 values
    n = 1;
    S(1) = N-I0; I(1) = I0; R(1) = 0;
    
    % -------------------- gillespie loop --------------------
    % continue loop if still in time frame and # infected is positive
    while t < tF-dt_avg && I(n) > 0
        
        % average number of newly infected per time unit
        w1 = beta*S(n)*I(n)/N;
        % average number of newly recovered at time, t
        w2 = gamma*I(n);
        % total number of new group changes at time, t
        W = w1 + w2;
        
        % Notes pertaining to the computation of dt: 
        %  * dt = -log(1-rand)/W <=> rand = 1-exp(-W*dt)
        %  * dt should be small enough so that movement bewteen states in
        %    unlikely to occur. What is the probability that there is no
        %    movement at any given time?
        %  * What is the probably that no one moved to I?
        %        P(not_I) = exp(-w1/W*t)
        %  * What is the probably that no one moved to R?
        %        P(not_I) = exp(-w2/W*t)
        %  * What is the probably that no one moved?
        %        P(not_I & not_R) = exp(-w1*t)*exp(-w2*t)
        %                         = exp(-(w1+w2)*t)
        %                         = exp(-W*t)
        %  * What is the probably that one moved?
        %        P(I or R) = 1 - P(not_I & not_R)
        %                  = 1 - exp(-W*t)
        %  * Probability Distribution
        %        P(t) *
        %             |*
        %   uniform   |   *            P(dt) = exp(-w*dt)
        %    random   |        *
        %    number   |               *
        %           --|----------------------*------> dt
        %
        %         solve for dt to get:   dt = -log(1-rand)/W
        %
        %  * dt is the average time to wait until the first move between
        %    states happens.
        
        % compute random time step, dt
        dt = -log(1-rand)/W;
        t = t + dt;
        
        time(n+1) = t;
        if rand < w1/W  % w1/W = probability of becoming infected at time, t
            S(n+1) = S(n) - 1;    % one leaves susceptible group
            I(n+1) = I(n) + 1;    % one enters infected group
            R(n+1) = R(n);        % recovered group unchanged
        else
            S(n+1) = S(n);        % susceptible group unchanged
            I(n+1) = I(n) - 1;    % one leaves infected group
            R(n+1) = R(n) + 1;    % one enters recovered group
        end
        n = n + 1;
    end
    
    % remove trailing zeros
    time(n+1:end) = [];
    S(n+1:end) = [];
    I(n+1:end) = [];
    R(n+1:end) = [];
    
    % update averages
    for k = 1:length(time)
        tk = time(k);
        [~,idx] = min(abs(time_std - tk));
        c(idx) = c(idx) + 1;
        S_avg(idx) = ((c(idx)-1)*S_avg(idx) + S(k))/c(idx);
        I_avg(idx) = ((c(idx)-1)*I_avg(idx) + I(k))/c(idx);
        R_avg(idx) = ((c(idx)-1)*R_avg(idx) + R(k))/c(idx);
    end
    
    % plot current averages
    pause(0.1)
    hold off;
    plot(time_std(S_avg ~= 0),S_avg(S_avg ~= 0),'x');
    hold on;
    plot(time_std(I_avg ~= 0),I_avg(I_avg ~= 0),'o');
    plot(time_std(R_avg ~= 0),R_avg(R_avg ~= 0),'+');
    
    % plot ode solution
    plot(t_ode,S_ode','-',t_ode,I_ode','-',t_ode,R_ode','-')
    
end
end

function dydt = odefcn(y,gamma,beta)

S = y(1); I = y(2); R = y(3);
N = S+I+R;
dSdt = -beta*S*I/N;
dIdt =  beta*S*I/N - gamma*I;
dRdt =  gamma*I;

dydt = [dSdt; dIdt; dRdt];
end