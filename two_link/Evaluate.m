function output = Evaluate(Index)

Handler_dynamics_HCg_model = SRD_get(['Handler_dynamics_HCg_model', num2str(Index)]);
Handler_t_linearization = SRD_get('Handler_reduced_dynamics_and_transverse_linearization');


p = [0 1]; % vrtl cnstr parameter vector
vrtl_cnstr = SRDt_VirtualConstraint(Handler_t_linearization.N_dof, length(p)-1, 'ordinary', ...
    Handler_t_linearization.c0, ...
    Handler_t_linearization.H0);


nmnl_trj = SRDt_get_nominal_trajectory(...
    'Handler_reduced_dynamics_and_transverse_linearization', Handler_t_linearization, ...    
    'p', [0 1], ...
    'vrtl_cnstr_obj', vrtl_cnstr, ...
    's0', 0.1, ...
    'dt', 0.01);


% Script for designing the controller

ReSolve = false;
if ReSolve
    
    N  = size(nmnl_trj.A, 3);
    dt = nmnl_trj.T/(N-1); % to consider zero
    tt = 0:dt:nmnl_trj.T;
    Qi = eye(3); Qi(1,1) = 1;
    Q = repmat(Qi,[1,1,N]);
    R = repmat(1,[1,1,N]);
    
    
    optns.solver    = 'sdpt3';
    optns.M         = 20;
    optns.sgma      = 0.1;
    optns.d         = 1e+5;
    optns.pol_type  = 'T';
    
    [K_mtrx, X] = solvePRDE(nmnl_trj.A,nmnl_trj.B,Q,R,tt,optns);
    
else
    K_mtrx = load('matlab.mat', 'K_mtrx');
end
% Nonlinear closed-loop

% options
no_iter = 5000;
delta_t = 1e-2;
add_disturbance = 1; % state distrubance at initial time

if add_disturbance
    delta_x0 = [0.01, 0.01, 0.01, 0.01];
else
    delta_x0 = zeros(1, 4);
end

% simulation time

x0_str = [nmnl_trj.q(:,1)' nmnl_trj.qd(:,1)'];
x0_dstbd = x0_str + delta_x0; % initial state

% integrator options
optns_ode45 = odeset('AbsTol',1e-10,'RelTol',1e-8);

% create ode fucntion to be used with ode45
% odefun = @(t,x,u) full_dnmcs_ode(x,u);
odefun = get_handler_full_dnmcs_ode(Handler_dynamics_HCg_model);
h_Intg2 = SRDt_get_handler_Intg2(p, vrtl_cnstr, Handler_t_linearization);

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
%     disp(['calculating ', num2str(i), ' out of ', num2str(no_iter)]);
    
    if rem( i, floor(no_iter / 100)) == 0
        disp(['calculating ', num2str( floor(100 * i / no_iter) ), ' %']);
    end
    
    t_span = [(i-1)*delta_t, i*delta_t];
    
    % Find current Poincare section;
    delta_x = x_dstbd(i,:) - [nmnl_trj.q' nmnl_trj.qd'];
    delta_x_nrm = vecnorm(delta_x, 2, 2); % find distance to each member of the set (nominal trajectory)
    [~, idx] = min(delta_x_nrm); % choose minimum distance
    
    % introduce the motion generator for readability of the code
    si = x_dstbd(i,1);
    sdi = x_dstbd(i,3);
    
       
    phi2i = Handler_t_linearization.H0*vrtl_cnstr.Phi(si, p);
    phi2_prmi = Handler_t_linearization.H0*vrtl_cnstr.Phi_prm(si, p);
    phi2_2prmi = Handler_t_linearization.H0*vrtl_cnstr.Phi_2prm(si, p);

    % Compute Transverse coordinates
    I_cur = h_Intg2(si, sdi, nmnl_trj.q(1, idx), nmnl_trj.qd(1, idx));
    y_cur = x_dstbd(i,2) - phi2i;
    y_d_cur = x_dstbd(i,4) - phi2_prmi*sdi;
    x_trsv_cur = [I_cur, y_cur, y_d_cur]';
   
    
    
    y(i) = y_cur;
    y_d(i) = y_d_cur;
    I(i) = I_cur;
    
    % Input that consistes of feedfowfard term and feedback terms
    u_ffrd = Handler_t_linearization.get_Uff(si, sdi, y_cur, y_d_cur, phi2i, phi2_prmi, phi2_2prmi);   
    u_fbck = K_mtrx(:,:,idx)*x_trsv_cur;
    u_cur = u_ffrd + 1/Handler_t_linearization.get_Nff(si, y_cur, phi2i, phi2_prmi) * u_fbck;
    
    u(i) = u_cur;
    ufb(i) = u_cur - u_ffrd;
    uff(i) = u_ffrd;
    
    % Integration of dynamics from t_{i} to t_{i+1} with u = u_{i}
    [~, x_cur] = ode45(@(t,x) odefun(t,x,u_cur),...
                       t_span, x_dstbd(i,:)', optns_ode45);                   
   
    x_dstbd(i+1,:) = x_cur(end,:);
end

output.x_cur = x_cur;
output.x_dstbd = x_dstbd;
output.y = y;
output.y_d = y_d;
output.I = I;

end







