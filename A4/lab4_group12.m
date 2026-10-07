
addpath(genpath(fullfile(pwd,'yalmip')));


S = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);

% clear local data
yalmip('clear');

% x = sdpvar(2,1);

c1 = 1/60;
c2 = 0.2;
c3 = 0.1;


u_vec = [0, 0.3, 0.8];  % Lockdown, distancing, no distancing

% set P = x1*x3?


A = [
    0 0 0 c1;
    0 -c2 0 0;
    0 c2 0 0;
    0 0 0 -c1;
];

B = [
    
];






