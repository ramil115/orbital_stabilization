close all; clear classes; clc;

g = SRD_get('Handler_dynamics_HCg_model');
h = SRD_get('Handler_reduced_dynamics_and_transverse_linearization');


p = [0 1]; % vrtl cnstr parameter vector
vrtl_cnstr = SRDt_VirtualConstraint(h.N_dof, length(p)-1, 'ordinary', h.c0, h.H0);


nmnl_trj = SRDt_get_nominal_trajectory(...
    'Handler_reduced_dynamics_and_transverse_linearization', h, ...    
    'p', [0 1], ...
    'vrtl_cnstr_obj', vrtl_cnstr, ...
    's0', 0.1, ...
    'dt', 0.01);


%% Script for designing the controller 
% N = size(nmnl_trj.A,3);
% dt = nmnl_trj.T/(N-1); % to consider zero 
% tt = 0:dt:nmnl_trj.T;
% Qi = eye(3); Qi(1,1) = 1;
% Q = repmat(Qi,[1,1,N]);
% R = repmat(1,[1,1,N]);
% 
% 
% optns.solver    = 'sdpt3';
% optns.M         = 20;
% optns.sgma      = 0.1; 
% optns.d         = 1e+5;
% optns.pol_type  = 'T';
% 
% [K_mtrx, X] = solvePRDE(nmnl_trj.A,nmnl_trj.B,Q,R,tt,optns);

load('matlab.mat')

%% Nonlinear closed-loop

% options
no_iter = 5000;
delta_t = 1e-2;
add_disturbance = 1; % state distrubance at initial time

if add_disturbance
    delta_x0 = [0.01, 0.01, 0.01, 0.01];
else
    delta_x0 = zeros(1,4);
end

% simulation time
t = 0:delta_t:no_iter*delta_t;

% Check feedforward consistency

uff_1 = zeros(1, length(nmnl_trj.s));
uff_2 = zeros(1, length(nmnl_trj.s));
for i = 1:length(nmnl_trj.s)
    uff_1(i) = pinv(g.get_T())*(g.get_H(nmnl_trj.q(:,i))*nmnl_trj.qdd(:,i) + g.get_C(nmnl_trj.q(:,i), nmnl_trj.qd(:,i))*nmnl_trj.qd(:,i) + ...
        g.get_g(nmnl_trj.q(:,i)));
    uff_2(i) = h.get_Uff(nmnl_trj.q(1,i), nmnl_trj.qd(1,i), 0, 0, h.H0*vrtl_cnstr.Phi(nmnl_trj.q(1,i), p),...
        h.H0*vrtl_cnstr.Phi_prm(nmnl_trj.q(1,i), p), h.H0*vrtl_cnstr.Phi_2prm(nmnl_trj.q(1,i), p));
end

if max(abs(uff_1-uff_2))< 1e-6
    disp('Feedforward is consistent')
else
    disp('Feedforward is NOT consistent')
end

x0_str = [nmnl_trj.q(:,1)' nmnl_trj.qd(:,1)'];
x0_dstbd = x0_str + delta_x0; % initial state

% integrator options
optns_ode45 = odeset('AbsTol',1e-10,'RelTol',1e-8);

% create ode fucntion to be used with ode45
% odefun = @(t,x,u) full_dnmcs_ode(x,u);
odefun = get_handler_full_dnmcs_ode(g);
h_Intg2 = SRDt_get_handler_Intg2(p, vrtl_cnstr, h);

% prelocate variables
x_dstbd = zeros(no_iter,4);
y = zeros(no_iter,1);
y_d = zeros(no_iter,1);
I = zeros(no_iter,1);
u = zeros(no_iter,1);
ufb = zeros(no_iter,1);
uff = zeros(no_iter,1);

x_dstbd(1,:) = x0_dstbd;

for i = 1:no_iter
    t_span = [(i-1)*delta_t, i*delta_t];
    
    % Find current Poincare section;
    delta_x = x_dstbd(i,:) - [nmnl_trj.q' nmnl_trj.qd'];
    delta_x_nrm = vecnorm(delta_x, 2, 2); % find distance to each member of the set (nominal trajectory)
    [~, idx] = min(delta_x_nrm); % choose minimum distance
    
    % introduce the motion generator for readability of the code
    si = x_dstbd(i,1);
    sdi = x_dstbd(i,3);
    
       
    phi2i = h.H0*vrtl_cnstr.Phi(si, p);
    phi2_prmi = h.H0*vrtl_cnstr.Phi_prm(si, p);
    phi2_2prmi = h.H0*vrtl_cnstr.Phi_2prm(si, p);

    % Compute Transverse coordinates
    I_cur = h_Intg2(si, sdi, nmnl_trj.q(1, idx), nmnl_trj.qd(1, idx));
    y_cur = x_dstbd(i,2) - phi2i;
    y_d_cur = x_dstbd(i,4) - phi2_prmi*sdi;
    x_trsv_cur = [I_cur, y_cur, y_d_cur]';
   
    
    
    y(i) = y_cur;
    y_d(i) = y_d_cur;
    I(i) = I_cur;
    
    % Input that consistes of feedfowfard term and feedback terms
    u_ffrd = h.get_Uff(si, sdi, y_cur, y_d_cur, phi2i, phi2_prmi, phi2_2prmi);   
    u_fbck = K_mtrx(:,:,idx)*x_trsv_cur;
    u_cur = u_ffrd + 1/h.get_Nff(si, y_cur, phi2i, phi2_prmi) * u_fbck;
    
    u(i) = u_cur;
    ufb(i) = u_cur - u_ffrd;
    uff(i) = u_ffrd;
    
    % Integration of dynamics from t_{i} to t_{i+1} with u = u_{i}
    [~, x_cur] = ode45(@(t,x) odefun(t,x,u_cur),...
                       t_span, x_dstbd(i,:)', optns_ode45);                   
   
    x_dstbd(i+1,:) = x_cur(end,:);
end

function h = get_handler_full_dnmcs_ode(g)
h = @(t,x,u) full_dnmcs_ode(x,u,g);
% h = @full_dnmcs_ode;
% function  dxdt = full_dnmcs_ode(x,u)
function  dxdt = full_dnmcs_ode(x,u,g)
    H = g.get_H(x(1:2));
    C = g.get_C(x(1:2), x(3:4));
    gr = g.get_g(x(1:2));        
    T = g.get_T();   
    dxdt = [x(3:4); H\(T * u - C * x(3:4) - gr)];
end
end













