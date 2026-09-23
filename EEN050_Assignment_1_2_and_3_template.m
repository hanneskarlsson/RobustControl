%% Robust and Nonlinear Control, EEN050:
% Template for assignment 1, 2, and 3
%-------------------------------------------------
%-------------------------------------------------
%Initialization
clear all;
close all;
clc
% READ THIS:
% - Run this section first and do not overwrite any of the variables here
% besides (if necessary) input names and output names for the ss-objects.

%Low fidelity F16 longitudinal model from Aircraft Control and Simulation (B.L.Stevens-F.L.Lewis) pp. 156
%The linearized dynamic of the airplane F16
A_n=[-0.127 -235 -32.2 -9.51 0.314;-7E-4 -0.969 0 0.908 -2E-4;0 0 0 1 0; 9E-4 -4.56 0 -1.58 0;0 0 0 0 -5];
%States
%[V(ft/s) speed, alpha(rad) angle of attack, theta(rad) pitch angle, q(rad/s) pitch rate, T(lb) engine_power]'
B_n=[0 -0.244;0 -0.00209; 0 0;10 -0.199; 1087 0];
%Control inputs
%[thrust (N); elevator_deflection(rad)]'
C_n=[0 57.3 0 0 0;0 0 0 1 0;0.0208 15.2 0 1.45 0];
%Measured outputs
%[alpha(deg); q(rad/s); normal_acceleration(ft/s^2) ]
D_n=[0 0;0 0;0 0.033];
%Note, elevator deflection has a direct effect on vertical/normal
%accelleration
Gn=ss(A_n,B_n,C_n,D_n);
tf(Gn)
Gn.InputName = 'utilde';
Gn.OutputName = 'y';
% Actuator dynamics (nominal):
GT = tf(1,[1/(2.5*10) (1/2.5+1/10) 1]); % Thrust
Ge = tf(1,[1/25 1]);                    % Elevator deflection
Ga  = ss([GT, 0;0, Ge]);                % Actuator dynamics
Ga.InputName = 'u';
Ga.OutputName = 'ytilde';
tf(Ga)

% Disturbance model:
dryden =  tf([0.9751 0.2491],[1 0.885 0.1958]);
Wd = ss([0;dryden]);
Wd.InputName = 'd';
Wd.OutputName = 'Wd';

% Noise filter Wn:
wn1 = rad2deg(0.001);
wn2 = 0.001;
ft2m = 0.3048; % feet/meter
wn3 = 0.001/ft2m;
Wn = ss(diag([wn1,wn2,wn3]));
Wn.OutputName = 'Wn';
Wn.InputName = 'n';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LQG design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the physical system
G = Gn*[Ga Wd];

[A,B,C,~] = ssdata(G); % Retrive model matrices
% Number of states, control inputs, process noise, measurement noise
nx = length(A); nu = 2; nw = 1; ny = 3;

Glqi = ss(A,B(:,1:2),C(1,:),0); % ss-object is used for the LQI design

% LQI design - weigthing matrices
Q = blkdiag(eye(nx),10000);
R = eye(nu);

% Compute the feedback gain K
[K,~,~] = lqi(Glqi,Q,R);

% Kalman filter design w. 'kalman'
Gkalman = ss(A,B(:,3),C,0);% ss-object for filter design
QN = eye(nw);  % Process noise covariance
RN = eye(ny);  % Measurement noise covariance
[~,L,~] = kalman(Gkalman,QN,RN); % Compute the filter

% Split the gain into the state feedback gain and the integral gain
Kx = K(:,1:end-1);  % state feedback
Ki = K(:,end);      % integral gain
Bu = B(:,1:2);      % this is the part of B that multiplies with u

% Construct the system matrices
Ac = [A-L*C-Bu*Kx -Bu*Ki;-C(1,:) 0];
Bc = [L zeros(nx,1);0 0 0 1];
Cc = -K;
Dc = 0;

% Represent the LQG controller as a dynamical system
% Its input is [ytilde;r]
% Its output is [u]
LQG = ss(Ac,Bc,Cc,Dc);
LQG.InputName={'ytilde(1)','ytilde(2)','ytilde(3)','r'};
LQG.Outputname = {'u(1)','u(2)'};
'You can ignore the above name conflict'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% A1/Ex1 - Open loop analysis
clc
G_F16 = Gn*Ga;                            % Airplane model
% Is the airplane model G_F16 is stable?
eigenvalues = eig(G_F16)

% all negative so it is stable

% Is the airplane model G_F16 minimal order?
[A,B,C,D] = ssdata(G_F16);
n = size(A,1);
rankc = rank(ctrb(A,B));
ranko = rank(obsv(A,C));
if n == rankc && n == ranko
    disp('System is minimal')
else
    disp('System is not minmal')
end
 % It is minimal


% Plot the singular values of the airplane model G_F16 (see the command
% 'sigma')
sigma(G_F16)

% 42.3 dB at 0.154 hz

% Compute the H-infinity norm of the airplane model
% (This can be done with the command 'hinfnorm'.)
G_F16Inf = mag2db(hinfnorm(G_F16))

% 42.3083 db

% Compute the H-2 norm of the airplane model 
% (The easiest way to do this is with the command 'norm'.)
% example: norm(ss-object of interest,2)

G_F16norm2 = norm(G_F16,2)
% Compute the H-2 norm of 'Gn' 
Gnnorm2 = norm(Gn,2)
pole(Gn)
tf(Gn)
% It is not strictly proper which causes it to have dirac delta a t=0 which
% has a infinite energy.
% 
%% A1/Ex2 - Uncertain actuator dynamics
% MATLAB-functions to use, 'tf', 'ureal', 'usample', 'step', 'sigma',
% 'blkdiag'.
clc
% Define the uncertainty parameters with 'ureal'
t1   =  ureal('t1',2.5,'percentage',30);
t2   =  ureal('t2',10.0,'percentage', 40);
t3   =  ureal('t3', 25,'percentage', 25);
gT   =  ureal('gT', 1.0, 'percentage',  50);
ge   =  ureal('ge', 1.0, 'percentage',  10);


% Define the uncertain dynamics for thrust
GTu=tf([gT],[1/(t1*t2),(1/t1+1/t2),1]);                % Thrust

% Define the uncertain dynamics for elevator deflection
Geu=tf(ge, [1/t3, 1])                               % Elevator deflection

% Define the uncertain actuator dynamics (block diagonal)
Gau = blkdiag(GTu,Geu);
Gau.InputName = 'u';
Gau.OutputName = 'ydelta';

% Create 100 samples of Gau1 and Gau2 with 'usample'
samples = 100; % number of samples
% it takes time to generate 100 samples so feel free to change this
% the value of 'samples' when you are figuring things out.

Geu_samples = usample(Geu,100); 
GTu_samples = usample(GTu,100);
%f
% Plot the step responses and bode diagrams (sigma plot)
% Feel free to use the following plotting routine:
figure(1)
subplot(2,2,1)
step(GTu_samples)
hold on
step(GTu.NominalValue,'red')
title('Thrust: step response')
legend('random samples','nominal')
 
subplot(2,2,2)
step(Geu_samples)
hold on
step(Geu.NominalValue,'red')
title('Elevator deflection: step response')
legend('random samples','nominal')
 
subplot(2,2,3)
sigma(GTu_samples)
hold on
sigma(GTu.NominalValue, 'red')
title('Thrust: singular values')
legend('random samples','nominal')

subplot(2,2,4)
sigma(GTu_samples)
hold on
sigma(GTu.NominalValue, 'red')
title('Elevator deflection: singular values')
legend('random samples','nominal')
hold off

% Both systems remain stable with the uncertaintiy
% The thursts magnitude is however affected
%
%% A1/Ex3 - (1)
% Simulation of LQG controller w. nominal actuator dynamics
% Consult the documentation for 'lsim' and 'connect' to better understand
% the code
clc
% Simulation parameters:
N = 1000;                  % number of time steps in simulation
T = linspace(0,50,N);      % time vector
flag_noise = 1;            % set to zero to remove noise
flag_x0 = 1;               % set to zero to initialize system at the origin

% Inputs:
r = zeros(N,1);            % reference signal
r(1:200)=-0.5; r(201:400)=1; r(401:600)=5; r(601:800)=-2; r(801:end)=0;
noise = randn(N,4);        % disturbance and measurement noise
U = [r noise*flag_noise];  % input vector

% Define the closed loop system using 'connect':
% Provide appropriate input/output names for each ss-object
Wd.InputName = 'd';
Wd.OutputName = 'Wd';

Ga.InputName = 'u';
Ga.OutputName = 'ydelta';
Wn.InputName = 'n';
Wn.OutputName = 'Wn';

Gn.InputName = 'utilde';
Gn.OutputName = 'y';

LQG.InputName = {'ytilde(1)','ytilde(2)','ytilde(3)','r'};
LQG.OutputName = 'u';

% Summation blocks:
Sum1 = sumblk('utilde = ydelta+Wd',2);
Sum2 = sumblk('ytilde = y+Wn',3);

% Choose inputs and outputs for the resulting closed loop system
inputs = {'r','n','d'};
outputs = {'y(1)'};     % y(1) corresponds with angle of attack [deg]

% Define the closed loop system:
LQG_clp = connect(Ga,Gn,Wd,Wn,LQG,Sum1,Sum2,inputs,outputs);

% Number of states in the closed loop system
nx_LQG = length(LQG_clp.A);
x0 = randn(nx_LQG,1)*flag_x0;   % initial state

% Simulating and plotting
Y_LQG = lsim(LQG_clp,U,T,x0);
figure(2)
plot(T,Y_LQG,T,r,'r--')
legend('angle of attack','reference signal')
title('Reference tracking w. LQG on the angle of attack')
ylim([-10, 10])

%% A1/Ex3 - (2)
% Closed loop simulation w. LQG and uncertain actuator dynamics
% Repeat the above simulation with the ss-object 'Ga' replaced by a sample
% of the uncertain ss-object for the actuator dynamics 'Gau'(c.f A1/Ex2).
% The input and output names of the samples actuator dynamics should
% correspond with the input and output names of 'Ga'.

% Take a random sample of the uncertain actuator dynamics
Ga_random = usample(Gau);
Ga_random.InputName = 'u';
Ga_random.OutputName = 'ydelta';

% Close the loop with Ga_random instead of Ga
clp_LQGu = connect(Ga_random,Gn,LQG,Wn,Wd,Sum1,Sum2,inputs,outputs);
% simulate and plot
Y_LQGu = lsim(clp_LQGu,U,T,x0);
figure(3)
plot(T,Y_LQGu,T,r,'r--')
legend('angle of attack','reference signal')
title('Reference tracking w. LQG on the angle of attack (random actuator dynamics)')
ylim([-10, 10])



%% A2/Ex1
clc
% Notes:
% Be mindful of the input/output dimensions when defining the filters.
% Feel free to use 'makeweight' when apropriate
% 
% Wn and Wd are defined in the first code section of this file.
% (Be aware of their input and output names!)
% To define the filter Wm (consult the documentation for 'ucover').
%
% Use 'connect' to contstruct the 'P' matrix.
% For convenience, order the inputs/outputs as follows
% Input: [udelta;r;n;d;u]
% Output: [ydelta;z;v] (v=[y+n;r])
%
% To compute the Hinf controller use 'hinfsyn'
% To compute the H2 controller use 'h2syn'

% Define the filters Wra, We, Wp, Wu


% Define WmT, Wme, and Wm (read the documentation for 'ucover')
% Feel free to use the samples taken in A1/Ex2
Wralphatf = tf(6.25^2,[1, 2*6.25, 6.25^2]);
Wra = ss([Wralphatf; 0; 0]);
Wra.OutputName = 'Wralpha'
Wra.InputName = 'r'

% We 
Wetf = calculateWeight(400, 4.3, 0.4);
We = ss([Wetf, 0, 0]);
We.OutputName = 'ze'
We.InputName = 'sum4'
%Wp
Wangle = calculateWeight(2.5,0.45, 0.015);
Wacc = calculateWeight(2.5,0.7,0.0063);
Wp = ss([Wangle,0 ,0 ;0 , 0, Wacc]);
Wp.OutputName = 'zp'
Wp.InputName = 'y'
%Wu
Wutf = tf(1/deg2rad(35)); % maybe this should not be rad since the w1 is in degree
Wu = ss([0,0 ; 0,Wutf]);
Wu.OutputName = 'zu'
Wu.InputName = 'utilde'


w = logspace(-2,2,500);

Geu_samples = usample(Geu,100); 
GTu_samples = usample(GTu,100);

Geu_samples_frd = frd(Geu_samples,w);
GTu_samples_frd = frd(GTu_samples,w);

[Wme, info1] = ucover(Geu_samples_frd, Ge,2, 'Outputmult');
[WmT, info2] = ucover(GTu_samples_frd, GT, 2, 'Outputmult');

%Wm = ss(blkdiag(WmT,Wme));
Wm = ss(blkdiag(info2.W1, info1.W1))
Wm.InputName  = 'udelta';
Wm.OutputName = 'Wm';

% Provide appropriate input/output names
Sum4 = sumblk('sum4 = Wralpha-y',3);
Sum5 = sumblk('utilde = ydelta+Wm+Wd',2);
Sum6 = sumblk('ytilde = Wn+y',3);
% Define the summation blocks (there are three):

% Define the appropriate inputs and outputs (the order matters!)
inputs = {'udelta','r','n','d','u'}
outputs = {'ydelta','ze','zp','zu','ytilde', 'r'}

%

P = connect(Ga, Gn, Wm, Wd, Wn, Wra, We, Wp, Wu, Sum4, Sum5, Sum6, inputs, outputs);
    
function W = calculateWeight(DCgain, wc, HFgain)
% calculateWeight Calculates K, z, and p for
%
%        W(s) = K * (s + z) / (s + p)
%
% Inputs:
%   DCgain - DC gain of the weighting function
%   wc     - crossover frequency [rad/s]
%   HFgain - high-frequency gain
%
% Outputs:
%   K - high-frequency gain
%   z - zero parameter
%   p - pole parameter

    % High-frequency gain
    K = HFgain;

    % From DC gain:
    % K*z/p = DCgain
    ratio = DCgain / K;

    % From crossover condition:
    % K^2*(wc^2 + z^2)/(wc^2 + p^2) = 1
    %
    % Since z = ratio*p:
    p = sqrt((wc^2*(1-K^2)) / (K^2*ratio^2 - 1));

    % Calculate z
    z = ratio*p;
    s = tf('s');
    W = K*((s + z)/(s + p));
end
%% A2/Ex2
% Compute the H-infinity controller

[Kinf, Ninf, gamma, info_inf] = hinfsyn(P, 4, 2);

% plot the singular values
figure(4)
sigma(Kinf)
grid on
title('H_inf: singular values')


% plot the singular values
figure(4)

%% A2/Ex3
% Compute the H2-controller
% Compute the H2-controller
[K_2, N_2, gamma_2, info_2] = h2syn(P, 4, 2);
% plot the singular values
figure(5)
sigma(K_2);
grid on
title('H_2: singular values');

% plot the singular values
figure(5)



%% A3/Ex1
% Nominal stability true
isstable(tf(Ninf))

% Nominal performance

%% ======================================================================
%  EEN050 - Assignment 3 solution
%  Paste these sections in place of the A3 skeleton in
%  EEN050_Assignment_1_2_and_3_template.m
%  (Assumes A1 and A2 sections have already been run, so that
%   Gn, Ga, Wd, Wn, Gau, LQG, P, Kinf, Ninf, gamma, info_inf, K_2, N_2
%   are already defined in the workspace.)
%  ======================================================================

%% A3/Ex1 - Closed loop analysis (H_infinity controller)
clc
w = logspace(-2,3,500);              % omega in [0.01, 1000] rad/s

N = Ninf;                            % N = lft(P,Kinf) from hinfsyn

% ---- Nominal Stability (NS) --------------------------------------------
NS = isstable(N);
fprintf('Nominal stability (NS): %d\n', NS)

% ---- Channel partition of N --------------------------------------------
% N inputs  (order): udelta(2)  r(1)  n(3)  d(1)
% N outputs (order): ydelta(2)  ze(1) zp(2) zu(2)
in_udelta  = 1:2;   in_r = 3;   in_n = 4:6;   in_d = 7;
out_ydelta = 1:2;   out_ze = 3; out_zp = 4:5; out_zu = 6:7;

N11 = N(out_ydelta, in_udelta);                          % udelta -> ydelta (RS)
N22 = N([out_ze out_zp out_zu], [in_r in_n in_d]);        % [r,n,d] -> [ze,zp,zu] (NP)

% ---- Nominal Performance (NP): sigma_bar(N22) < 1 -----------------------
figure
sigma(N22, w); grid on
title('NP test:  \sigma(N_{22})  (must stay below 0 dB)')

% ---- Robust Stability (RS): sigma_bar(N11) < 1 ---------------------------
figure
sigma(N11, w); grid on
title('RS test:  \sigma(N_{11})  (must stay below 0 dB)')

% ---- Robust Performance (RP), quick sufficient test ----------------------
% sigma_bar(N) < 1 for all w  (equivalent to checking gamma<1 pointwise;
% conservative because it treats the performance block and Delta_m as one
% single full block instead of two separate blocks)
figure
sigma(N, w); grid on
title('RP test (sufficient):  \sigma(N)  (must stay below 0 dB)')

% ---- Robust Performance (RP), exact test via mu (structured singular value)
% Delta_hat = blkdiag(Delta_m , Delta_p)
%   Delta_m : 2x2 full complex block   (the actual uncertainty)
%   Delta_p : 5x5 full complex block   (fictitious performance block,
%             size = dim([ze;zp;zu]) = 1+2+2 = 5)
blk = [2 2; 5 5];
Nfrd = frd(N, w);
[mubnds, ~] = mussv(Nfrd, blk);

% mubnds is an frd object: mubnds(1,1) = upper bound, mubnds(1,2) = lower bound
muUpper = squeeze(mubnds(1,1).ResponseData);   % convert to a plain double vector
muLower = squeeze(mubnds(1,2).ResponseData);

figure
semilogx(w, 20*log10(muUpper)); hold on
semilogx(w, 20*log10(muLower));
yline(0,'r--')
grid on
xlabel('\omega [rad/s]'); ylabel('\mu bounds [dB]')
title('RP test (exact, via \mu)')
legend('\mu upper bound','\mu lower bound','0 dB')
% --- Report the results ---
disp('NS holds if isstable(N) = 1')
disp('NP holds if sigma_bar(N22) stays below 0 dB for all w in [0.01,1000]')
disp('RS holds if sigma_bar(N11) stays below 0 dB for all w in [0.01,1000]')
disp('RP holds if mu(N) (or, conservatively, sigma_bar(N)) stays below 0 dB')

%% A3/Ex2 - Closed loop simulation: LQG vs Hinf vs H2
clc

% Give Kinf and K_2 the SAME input/output names as LQG so they can be
% dropped into exactly the same feedback structure as in A1/Ex3
Kinf.InputName  = {'ytilde(1)','ytilde(2)','ytilde(3)','r'};
Kinf.OutputName = 'u';
K_2.InputName   = {'ytilde(1)','ytilde(2)','ytilde(3)','r'};
K_2.OutputName  = 'u';

% One random sample of the uncertain actuator dynamics.
% Use the SAME sample for all three controllers so the comparison is fair.
Ga_random = usample(Gau);
Ga_random.InputName  = 'u';
Ga_random.OutputName = 'ydelta';

% Re-affirm the names of the fixed blocks (as in A1/Ex3)
Gn.InputName = 'utilde';  Gn.OutputName = 'y';
Wd.InputName = 'd';       Wd.OutputName = 'Wd';
Wn.InputName = 'n';       Wn.OutputName = 'Wn';

Sum1 = sumblk('utilde = ydelta+Wd',2);
Sum2 = sumblk('ytilde = y+Wn',3);

inputs  = {'r','n','d'};
outputs = {'y(1)'};        % angle of attack [deg]

% Closed loop systems (same plant + Ga_random, three different controllers)
LQG_clp  = connect(Ga_random,Gn,Wd,Wn,LQG, Sum1,Sum2,inputs,outputs);
Hinf_clp = connect(Ga_random,Gn,Wd,Wn,Kinf,Sum1,Sum2,inputs,outputs);
H2_clp   = connect(Ga_random,Gn,Wd,Wn,K_2, Sum1,Sum2,inputs,outputs);

% Simulation parameters
N_sim = 1000;
T = linspace(0,50,N_sim);
flag_noise = 1;
flag_x0 = 1;

r = zeros(N_sim,1);
r(1:200)=-0.5; r(201:400)=1; r(401:600)=5; r(601:800)=-2; r(801:end)=0;
noise = randn(N_sim,4);
U = [r noise*flag_noise];

% Initial states (each closed loop has its own state dimension)
nx_LQG  = length(LQG_clp.A);
nx_Hinf = length(Hinf_clp.A);
nx_H2   = length(H2_clp.A);

x0_LQG  = randn(nx_LQG,1) *flag_x0;
x0_Hinf = randn(nx_Hinf,1)*flag_x0;
x0_H2   = randn(nx_H2,1)  *flag_x0;

% Simulate
YLQG = lsim(LQG_clp, U, T, x0_LQG);
Yinf = lsim(Hinf_clp,U, T, x0_Hinf);
Y2   = lsim(H2_clp,  U, T, x0_H2);

figure(8)
plot(T,Yinf,T,Y2,T,YLQG,T,r,'r--')
title('Reference tracking on angle of attack [deg] (uncertain actuator dynamics)')
xlabel('time [s]'); ylabel('\alpha [deg]')
legend('H_\infty','H_2','LQG','reference')
ylim([-3,6])
grid on

%% A3/Ex3 - Conclusion (discussion points for the report)
% Use the plot from A3/Ex2 (and, ideally, repeat it for a few different
% random samples of Ga to see if the ranking is consistent) to comment on:
%
% 1. Which controller tracks the reference fastest / with least overshoot?
% 2. Which controller is least sensitive to the randomly sampled actuator
%    uncertainty (compare against the nominal-Ga simulation from A1/Ex3)?
% 3. Which controller rejects noise (n) and disturbance (d) best?
% 4. Is it a fair comparison?
%    - LQG's Q,R weights were tuned "by feel" for the nominal plant only;
%      it was never told anything about the actuator uncertainty.
%    - Hinf/H2 were explicitly synthesized against Wm (which covers the
%      100 parametric samples), so they "know about" the uncertainty by
%      construction - an intrinsic advantage in this comparison.
%    - The three designs also don't share a common performance
%      specification (LQG uses quadratic state/ input cost with an
%      integrator; Hinf/H2 use frequency-weighted tracking error/ output/
%      input weights) so "better" partly reflects different design intent,
%      not just different mathematical techniques.
%    - Only ONE random Ga sample is used per run; a fair robustness
%      comparison would look at the worst case (or many samples) rather
%      than a single draw.

%% A3/Ex2 Simulation
% - Construct with 'connect' the  closed loop system between the LQG, Hinf,
%   and H2 controllers and the open loop plant.
% The open loop plant consists of:
% - Ga (nominal value or a sample of the uncertain actuator dynamics)
% - Wn, Wd, and, Gn
% - The delta block and the filters Wra, We, Wp, Wu, Wm are not included!

% Take a random sample of the uncertain actuator dynamics

% Define the appropriate summation blocks (there are two)


% Choose the apropriate inputs and outputs:


% Define the closed loop systems:
% LQG closed loop system
%LQG_clp = connect(...);

% H-infinity closed loop system
%Hinf_clp = connect(...);

% H-2 closed loop system
%H2_clp = connect(...);


% Feel free to use the following plotting routine:
% Simulation parameters:
% N = 1000;                  % number of time steps in simulation
% T = linspace(0,50,N);      % time vector
% flag_noise = 1;            % set to zero to remove noise
% flag_x0 = 1;               % set to zero to initialize system at the origin
% % Inputs:
% r = zeros(N,1);            % reference signal
% r(1:200)=-0.5; r(201:400)=1; r(401:600)=5; r(601:800)=-2; r(801:end)=0;
% noise = randn(N,4);        % Disturbance and measurement noise
% U = [r noise*flag_noise];  % Input vector
% 
% % Initial state(s)
% % Number of states is the LQG closed loop system
% nx_LQG = length(LQG_clp.A);
% % Number of states in the H_(2,infinity) closed loop systems
% nx_H = length(Hinf_clp.A);
% 
% x0_LQG = randn(nx_LQG,1)*flag_x0;   % initial state (LQG)
% x0_H = randn(nx_H,1)*flag_x0;       % initial state (Hinf,H2)
% 
% % Simulate:
% Yinf =  lsim(Hinf_clp,U,T,x0_H);
% Y2   =  lsim(H2_clp,U,T,x0_H);
% YLQG =  lsim(LQG_clp,U,T,x0_LQG);
% 
% % Feel free to use the following plotting routine.
% figure(8)
% plot(T,Yinf,T,Y2,T,YLQG,T,r,'r--')
% title('Reference tracking on angle of attack [deg]')
% legend('Hinf','H2','LQG','reference')
% ylim([-3,6])