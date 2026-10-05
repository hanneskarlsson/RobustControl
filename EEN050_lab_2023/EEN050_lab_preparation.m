warning off
format compact
clc
clear all
close all

colorOrder = get(gca,'colorOrder');
%% SETUP_SRV02_EXP14_2D_GANTRY
    %
    % Sets the necessary parameters to run the SRV02 Experiment #13: Position
    % Control of 2-DOF Robot laboratory using the "s_srv02_2d_robot" and
    % "q_srv02_2d_robot" Simulink diagrams.
    % 
    % Copyright (C) 2008 Quanser Consulting Inc.
    %
    clear all;
    qc_get_step_size =1/1000; %0.001;
    deltaT = qc_get_step_size;
    Ts = qc_get_step_size;
    nPoints = (50 * 3000) -1;
%
%% Initialization Settings(DONT CHANGE ANYTHING !!)
    EXT_GEAR_CONFIG = 'HIGH';
    ENCODER_TYPE = 'E';
    TACH_OPTION = 'YES';
    LOAD_TYPE = 'NONE';
    K_AMP = 1;
    AMP_TYPE = 'VoltPAQ';
    VMAX_DAC = 10;
    ROTPEN_OPTION = '2DGANTRY-E';
    PEND_TYPE = 'MEDIUM_12IN';
    THETA_MAX = 35 * pi/180;
    ALPHA_MAX = 15.0 * pi/180;
    CONTROL_TYPE = 'AUTO';   
    X0 = pi/180*[0, 0, 0, 0];
    [ Rm, kt, km, Kg, eta_g, Beq, Jm, Jeq_noload, eta_m, K_POT, K_TACH, K_ENC, VMAX_AMP, IMAX_AMP ] = config_srv02( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, AMP_TYPE, LOAD_TYPE );
    [ g, mp, Lp, lp, Jp_cm, Bp, RtpnOp, RtpnOff, K_POT_PEN ] = config_sp( PEND_TYPE, ROTPEN_OPTION );
    [ Lb, Jarm, K_POT_2DP, K_ENC_2DP ] = config_2d_gantry( Jeq_noload );
    K_ENC_2DIP = [-1,1].*K_ENC_2DP;
    wcf_1 = 2 * pi * 5;
    zetaf_1 = 0.9;
    wcf_2 = wcf_1;
    zetaf_2 = zetaf_1;

    %
%%%
%%%%
%%%%%
%%%%%% DO NOT CHANGE ANYTHING ABOVE THIS AREA !! Place your code below.
%%%%
%%%
%

%% Exercise 1

% Define the uncertain parameters Mp, Lp, Jp, and Co 
% using the command "ureal"


Mp = ureal('Mp', 0.1270, 'percentage', 50);
Lp = ureal('Lp', 0.3111, 'percentage',50);
Jp = ureal('Jp', 0.0012, 'percentage', 50);
Co = ureal('Co', 0.1285,'percentage', 50);



Lr = 0.1270; 
theta =  0;
alpha =  0;
dtheta =  0;
dalpha =  0;
Jr = 0.0083;
Dr = 0.0690;
g = 9.810;

Adelta = [0         0         1         0;
      0         0         0         1;
      0   ((Lr*Lp^3*Mp^2*dtheta^2)/2 + Lr*g*Lp^2*Mp^2)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - ((dalpha*dtheta*Lp^4*Mp^2)/2 + 2*Jp*dalpha*dtheta*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),- ((Jp*(4*Mp*alpha*dalpha*Lp^2 + 8*Dr))/2 + (Lp^4*Mp^2*alpha*dalpha)/2)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (- Lr*alpha*dtheta*Lp^3*Mp^2 + Dr*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),-((alpha*dtheta*Lp^4*Mp^2)/2 + 2*Jp*alpha*dtheta*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr));
      0   (dtheta^2*(Lp^2*Lr^2*Mp^2 + Jr*Lp^2*Mp) + 2*Lp*Lr^2*Mp^2*g + 2*Jr*Lp*Mp*g)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (Lp^3*Lr*Mp^2*dalpha*dtheta)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)), (2*alpha*dtheta*(Lp^2*Lr^2*Mp^2 + Jr*Lp^2*Mp))/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (Lr*alpha*dalpha*Lp^3*Mp^2 + 2*Dr*Lr*Lp*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),-(Lp^3*Lr*Mp^2*alpha*dtheta)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) ];
  
Bdelta = [0; 
      0;
      (4*Co*Jp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) + (Co*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr));
      (2*Co*Lp*Lr*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr))];

  Cdelta = [1  0  0  0;
     0  1  0  0];
 
 Ddelta = [0;0];


AdeltaNominal = Adelta.NominalValue
[vec, vals]= eig(AdeltaNominal)

%% Exercise 2
P_nom =  ss(Adelta.NominalValue,Bdelta.NominalValue, Cdelta, Ddelta);

Punc = ss(Adelta,Bdelta, Cdelta, Ddelta);
N = 100;
Psamples = usample(Punc, N);

w = logspace(-2,2,500);

Psample_frd = frd(Psamples,w);

[P, info] = ucover(Psample_frd, P_nom ,4, 'InputMult');
wim = info.W1;
WiM =10* blkdiag(wim, wim);


%% Exercise 3
Wu = ss(0.5 * eye(2));
wr = tf(1,[1,1])
Wr = ss(blkdiag(wr,wr));


tracking_err = [0.0025, 0.0251, 0.001];
wpVals = 1./(2*tracking_err)

wp1 = wpVals(1);
Wp = ss(blkdiag(wp1, wp1,wp1, wp1))

noiseAmpRad = 0.3*pi/180;
Wn = ss(noiseAmpRad * eye(4));

Pnom =ss(blkdiag(P_nom,P_nom));
Pnom.InputName = {'v1','v2'};
Pnom.OutputName = {'y1', 'y2', 'y3', 'y4'};

WiM = ss(WiM);
WiM.InputName = {'u1','u2'}
WiM.OutputName = {'yD1', 'yD2'}

Wu.InputName = {'u1', 'u2'}
Wu.OutputName  = {'zu1', 'zu2'}

Wr.InputName = {'r1', 'r2'}
Wr.OutputName = {'rf1', 'rf2'}

Wp.InputName = {'e1', 'e2', 'e3', 'e4'}
Wp.OutputName = {'ze1', 'ze2', 'ze3', 'ze4'}

Wn.InputName = {'n1', 'n2', 'n3', 'n4'}
Wn.OutputName  = {'wn1', 'wn2', 'wn3', 'wn4'}

Sv1 = sumblk('v1 = u1 + uD1');
Sv2 = sumblk('v2 = u2 + uD2');

Se1 = sumblk('e1 = y1 - rf1');  
Se2 = sumblk('e2 = y2');        
Se3 = sumblk('e3 = y3 - rf2');  
Se4 = sumblk('e4 = y4');        

Sm1 = sumblk('ym1 = y1 + wn1');
Sm2 = sumblk('ym2 = y2 + wn2');
Sm3 = sumblk('ym3 = y3 + wn3');
Sm4 = sumblk('ym4 = y4 + wn4');

PaugInputs = { ...
    'uD1','uD2', ...                  
    'r1','r2', ...                    
    'n1','n2','n3','n4', ...          
    'u1','u2'};                       

PaugOutputs = { ...
    'yD1','yD2', ...                  
    'zu1','zu2', ...                  
    'ze1','ze2','ze3','ze4', ...      
    'ym1','ym2','ym3','ym4'};         

Paug = connect( ...
    Pnom,WiM,Wu,Wr,Wp,Wn, ...
    Sv1,Sv2,Se1,Se2,Se3,Se4, ...
    Sm1,Sm2,Sm3,Sm4, ...
    PaugInputs,PaugOutputs);



%% Exercise 4
% Design and compute the controller, and call it Chinf
% Chinf = hinfsyn(...)
[Chinf,CLnom,gamma] = hinfsyn(Paug, 4,2);

fprintf('Achieved H-infinity gamma: %.6g\n',gamma);
fprintf('Nominal closed-loop stable: %d\n',isstable(CLnom));

w = logspace(-2,3,500);              % omega in [0.01, 1000] rad/s

  
% input udelta 1:2 r= 3 n = 4:6 d = 7
% output ydelta = 1:2 ze = 3 zp = 4:5 zu = 6:7 
% Nominal performance

figure(6)
N11 = CLnom(1:2, 1:2);
sigma(N11, w) % this should be less than 1 0db for 
yline(0,'r--')
title('RS test:  \sigma(N_{11})  (must stay below 0 dB)')
fprintf('Achieved H-infinity N11: %.6g\n',norm(N11, inf) );
   
figure (8)
sigma(CLnom, w)
yline(0,'r--')
title('RP test:  \sigma  (must stay below 0 dB)')
% also not stable sincce we dont have nominal performance





%% NB: To run the simulation, the nominal model has to be in the workspace with the variable name "Pnom"
%
%%%
%%%%
%%%%%
%%%%%% Closed loop simulation environment
%%%%
%%%
%
% NB: To run the simulation, the nominal model has to be loaded to the
% workspace with the variable name "Pnom".

[ah,bh,ch,dh] = ssdata(Chinf);


figure(1)
clf;
simTime = 10;
xinit=(pi/180)*[3 3 0 0 3 3 0 0];

try
sim('Simhinf.slx')
subplot(2,1,1)
hold on
    plot(simStates.Time,simStates.Data(:,1),'linewidth',2)
    plot(simStates.Time,simStates.Data(:,3),'linewidth',2)    
    plot(simStates.Time,simStates.Data(:,2),'--','linewidth',2)
    plot(simStates.Time,simStates.Data(:,4),'--','linewidth',2)
    legend('thetaX','alphaX','thetaY','alphaY','NthetaX','NalphaX','NthetaY','NalphaY')
    ylabel('Angle [DEG]')

subplot(2,1,2)
hold on
    plot(simVoltage.Time, simVoltage.Data(:,1),'linewidth',2)
    plot(simVoltage.Time, simVoltage.Data(:,2),'--','linewidth' ,2)
    legend('Voltage X','Voltage Y')
    ylabel('Voltage [V]')
    axis([0 simTime -11 11])
catch e
    disp('Simulation failed')
end
OCL=1;

%%
%load('labb_LQR_controller.mat')
%open('ExperimentRobustControl')
