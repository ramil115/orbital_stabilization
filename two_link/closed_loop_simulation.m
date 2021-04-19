% Define number of DoF
% s = c0*q // motion generator as a linear combination of q
% y = H0*q - phi(s) // define y in the following way
% p is parameter vector for polynom/bezier curve 
N_dof = 2;
c0 = [1, 0];
H0 = [0, 1];
p = [0 1];

vrtl_cnstr = VirtualConstraint(N_dof, length(p)-1,...
        'ordinary', c0, H0);


[s_str, sd_str, q_str, qd_str, qdd_str, T, A, B] = SRDt_get_nominal_trajectory(...
    'N_dof', 2, ...
    'c0', [1, 0], ...
    'H0', [0, 1], ...
    'p', [0 1], ...
    's0', 0.1, ...
    'dt', 0.01);
    
   



%% Script for designing the controller 
N = size(A,3);
dt = T/(N-1); % to consider zero 
tt = 0:dt:T;
Qi = eye(3); Qi(1,1) = 1;
Q = repmat(Qi,[1,1,N]);
R = repmat(1,[1,1,N]);


optns.solver    = 'sdpt3';
optns.M         = 20;
optns.sgma      = 0.1; 
optns.d         = 1e+5;
optns.pol_type  = 'T';

[K_mtrx, X] = solvePRDE(A,B,Q,R,tt,optns);


%% Nonlinear closed-loop

% options
no_iter = 5000;
delta_t = 1e-2;
add_disturbance = 1;

if add_disturbance
    delta_x0 = [0.01, 0.01, 0.01, 0.01];
else
    delta_x0 = zeros(1,4);
end

% simulation time
t = 0:delta_t:no_iter*delta_t;

% Check feedforward consistency

uff_1 = zeros(1, length(s_str));
uff_2 = zeros(1, length(s_str));
for i = 1:length(s_str)
    uff_1(i) = [1 0]*(g_dynamics_H(q_str(:,i))*qdd_str(:,i) + g_dynamics_C(q_str(:,i), qd_str(:,i))*qd_str(:,i) + ...
        g_dynamics_g(q_str(:,i)));
    uff_2(i) = get_U_ff(q_str(1,i), qd_str(1,i), 0, 0, H0*vrtl_cnstr.Phi(q_str(1,i), p),...
        H0*vrtl_cnstr.Phi_prm(q_str(1,i), p), H0*vrtl_cnstr.Phi_2prm(q_str(1,i), p));
end




x0_str = [q_str(:,1)' qd_str(:,1)'];
x0_dstbd = x0_str + delta_x0; % initial state

% integrator options
optns_ode45 = odeset('AbsTol',1e-10,'RelTol',1e-8);

% create ode fucntion to be used with ode45
odefun = @(t,x,u) full_dnmcs_ode(x,u);


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
    delta_x = x_dstbd(i,:) - [q_str' qd_str'];
    delta_x_nrm = vecnorm(delta_x, 2, 2); % find distance to each memeber of the set (nominal trajectory)
    [~, idx] = min(delta_x_nrm); % choose minimum distance
    
    % introduce the motion generator for readability of the code
    si = x_dstbd(i,1);
    sdi = x_dstbd(i,3);
    
       
    phi2i = H0*vrtl_cnstr.Phi(si, p);
    phi2_prmi = H0*vrtl_cnstr.Phi_prm(si, p);
    phi2_2prmi = H0*vrtl_cnstr.Phi_2prm(si, p);

    % Compute Transverse coordinates
    I_cur = Intg2(si, sdi, q_str(1, idx), qd_str(1, idx), p, vrtl_cnstr);
    y_cur = x_dstbd(i,2) - phi2i;
    y_d_cur = x_dstbd(i,4) - phi2_prmi*sdi;
    x_trsv_cur = [I_cur, y_cur, y_d_cur]';
   
    
    
    y(i) = y_cur;
    y_d(i) = y_d_cur;
    I(i) = I_cur;
    
    % Input that consistes of feedfowfard term and feedback terms
    u_ffrd = get_U_ff(si, sdi, y_cur, y_d_cur, phi2i, phi2_prmi, phi2_2prmi);   
    u_fbck = K_mtrx(:,:,idx)*x_trsv_cur;
    u_cur = u_ffrd + 1/get_N(si, y_cur, phi2i, phi2_prmi) * u_fbck;
    
    u(i) = u_cur;
    ufb(i) = u_cur - u_ffrd;
    uff(i) = u_ffrd;
    
    % Integration of dynamics from t_{i} to t_{i+1} with u = u_{i}
    [~, x_cur] = ode45(@(t,x) odefun(t,x,u_cur),...
                       t_span, x_dstbd(i,:)', optns_ode45);                   
   
    x_dstbd(i+1,:) = x_cur(end,:);
end


%%

function z_dot = reduced_dynamics_ode(t, z, Phi, Phi_prm, Phi_2prm)
     z_dot = [z(2); -1/get_alpha(Phi, Phi_prm)*...
                    (get_gamma(Phi) + ...
                        get_beta(Phi, Phi_prm, Phi_2prm)*z(2)^2)];
end


function  dxdt = full_dnmcs_ode(x,u)
    M = g_dynamics_H(x(1:2));
    C = g_dynamics_C(x(1:2), x(3:4));  
    G = g_dynamics_g(x(1:2));
        
    B = [1; 0];
   
    dxdt = [x(3:4); M\(B * u - C * x(3:4) - G)];
end

function  [s_str, sd_str, q_str, qd_str, qdd_str, T] = get_nominal_trajectory(varargin)
    Parser = inputParser;
    Parser.FunctionName = 'get_nominal_trajectory';
    Parser.addOptional('N_dof', []);
    Parser.addOptional('c0', []); % s = c0*q
    Parser.addOptional('H0', []); % H0*q = phi(s)
    Parser.addOptional('p', []); % virtual constraints parameter vector, e.g. phi(s) = p(1) + p(2)*s...
    Parser.addOptional('s0', []); % initial value of motion generator
    Parser.addOptional('dt', []); % discretization time
    
    Parser.parse(varargin{:});
    
    p = Parser.Results.p;

    vrtl_cnstr = VirtualConstraint(Parser.Results.N_dof, length(Parser.Results.p)-1,...
        'ordinary', Parser.Results.c0, Parser.Results.H0);
    
    %     func1 = @(s) get_alpha(vrtl_cnstr.Phi(s, p), vrtl_cnstr.Phi_prm(s, p));
    %     func2 = @(s) get_gamma(vrtl_cnstr.Phi(s, p));
    %     s_interval = [-10 10];
    % fplot(func1,s_interval)
    % roots1 = fzero(func1,0)
    % fplot(func2,s_interval)
    % roots2 = fzero(func2,0)

    tspan= 0:Parser.Results.dt:10;
    x0 = [Parser.Results.s0;0];
    optns = odeset('RelTol',1e-9,'AbsTol',1e-9,'NormControl','on');

    [t, x] = ode45( @(t,x)reduced_dynamics_ode(t,x,vrtl_cnstr.Phi(x(1), p), ...
        vrtl_cnstr.Phi_prm(x(1), p),...
        vrtl_cnstr.Phi_2prm(x(1), p)), tspan, x0, optns);

    [~,locs] = findpeaks(x(:,1));
    T = t(locs(1)); 
    T_ind = locs(1);    
    
    s_str = x(1:T_ind,1);
    sd_str = x(1:T_ind,2);
    
    sdd_str = zeros(T_ind,1);

    for i=1:length(T_ind)
        sdd_str(i) = -1/get_alpha(vrtl_cnstr.Phi(x(i,1), p), vrtl_cnstr.Phi_prm(x(i,1), p))*...
                  (get_gamma(vrtl_cnstr.Phi(x(i,1), p)) + ...
                  get_beta(vrtl_cnstr.Phi(x(i,1), p), vrtl_cnstr.Phi_prm(x(i,1), p), ...
                  vrtl_cnstr.Phi_2prm(x(i,1), p))*x(i,2)^2);
    end

q_str = zeros(2, T_ind);
qd_str = zeros(2, T_ind);
qdd_str = zeros(2, T_ind);
for i = 1:T_ind
    q_str(:,i) = vrtl_cnstr.Phi(s_str(i), p);
    qd_str(:,i) = vrtl_cnstr.Phi_prm(s_str(i), p)*sd_str(i);
    qdd_str(:,i) = vrtl_cnstr.Phi_2prm(s_str(i), p)*sd_str(i)^2 + ...
        vrtl_cnstr.Phi_prm(s_str(i), p)*sdd_str(i);
end

end

function  [A, B] = get_AB_matrices(varargin)
    Parser = inputParser;
    Parser.FunctionName = 'get_AB_matrices';
    Parser.addOptional('N_dof', []);
    Parser.addOptional('H0', []); % H0*q = phi(s)
    Parser.addOptional('c0', []); % s = c0*q
    Parser.addOptional('p', []); % virtual constraints parameter vector, e.g. phi(s) = p(1) + p(2)*s...
    Parser.addOptional('s_str', []); % initial value of motion generator
    Parser.addOptional('sd_str', []); % discretization time
    
    Parser.parse(varargin{:});
    
    vrtl_cnstr = VirtualConstraint(Parser.Results.N_dof, length(Parser.Results.p)-1,...
    'ordinary', Parser.Results.c0, Parser.Results.H0);
    
    p = Parser.Results.p;
    s_str = Parser.Results.s_str;
    sd_str = Parser.Results.sd_str;
    
    A = zeros(2*Parser.Results.N_dof-1,2*Parser.Results.N_dof-1,length(s_str));
    B = zeros(2*Parser.Results.N_dof-1,Parser.Results.N_dof-1,length(s_str));

    for i = 1:length(s_str)
        s_i = s_str(i);
        sd_i = sd_str(i);
        phi_i = Parser.Results.H0*vrtl_cnstr.Phi(s_i, p);
        phi_prm_i = Parser.Results.H0*vrtl_cnstr.Phi_prm(s_i, p);
        phi_2prm_i = Parser.Results.H0*vrtl_cnstr.Phi_2prm(s_i, p);
        A(:,:,i) = get_A_transv(s_i,sd_i,phi_i,phi_prm_i,phi_2prm_i);
        B(:,:,i) = get_B_transv(s_i,sd_i,phi_i,phi_2prm_i);
    end
    
end





