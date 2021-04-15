close all; clear classes; clc;

InitialPosition = SRD_get('InitialPosition');

Handler_State = SRD_get_handler__state('InitialPosition', InitialPosition, ...
    'InitialVelocity', zeros(size(InitialPosition)));

Handler_IK_Solution = SRD_get('Handler_IK_Solution');

Handler_Time = SRDHandler_Time();

Handler_Desired_State = SRD_get_handler__desired_state(...
    'Handler_ControlInput', Handler_IK_Solution, ...
    'Handler_Time',   Handler_Time);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model = SRD_get('Handler_dynamics_Linearized_Model');
Handler_Constraints_Model = SRD_get('Handler_Constraints_Model');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = Handler_IK_Solution.TimeExpiration;
% tf = 0.4;
time_table = 0:0.01:tf;

n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;


% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model


[A_table, B_table, c_table, x_table, u_table, dx_table] = ...
    SRD_LinearModel_GenerateTable('Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_dynamics_Linearized_Model', Handler_dynamics_Linearized_Model, ...
    'Handler_IK_Solution', Handler_IK_Solution, ...
    'TimeTable', time_table);

% n_constrained = 5;
[N_table, G_table, F_table] = SRD_ConstraintsModel_GenerateTable(...
    'Handler_Constraints_Model', Handler_Constraints_Model, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'x_table', x_table, ...
    'new_dimentions', []);
% n_constrained = size(N_table, 2);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controllers

Q = 100*eye(2 * n);
R = 0.01*eye(Handler_dynamics_generalized_coordinates_model.dof_control);
Count = size(A_table, 3);
K_table = SRD_LQR_GenerateTable(A_table, B_table, repmat(Q, [1, 1, Count]), repmat(R, [1, 1, Count]));


% K_table = SRD_CTC_GenerateTable('Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
%     'Handler_IK_Solution', Handler_IK_Solution, ...
%     'Kp', 200*eye(n), 'Kd', 100*eye(n), 'TimeTable', time_table);


% Q = 100*eye(2 * n);
% R = 0.01*eye(Handler_dynamics_generalized_coordinates_model.dof_control);
% Count = size(A_table, 3);
% K_table = SRD_CLQR_GenerateTable(A_table, B_table, repmat(Q, [1, 1, Count]), repmat(R, [1, 1, Count]), N_table);


% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA_table, cc_table] = SRD_CloseLoop_GenerateTable(A_table, B_table, c_table, K_table, x_table, u_table);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ode_fnc_handle = SRD_get_ode_fnc_from_ClosedLoopLinearSystem(AA_table, cc_table, time_table);

ode_fnc_handle = SRD_get_ode_fnc_from_ClosedLoopConstrainedLinearSystem...
    (AA_table, cc_table, G_table, F_table, time_table);


x0 = [InitialPosition; zeros(size(InitialPosition))];

[time_table_0, solution_tape] = ode45(ode_fnc_handle, [0, tf], x0);

figure('Color', 'w')
plot(time_table_0, solution_tape, 'LineWidth', 3); hold on;
plot(time_table', x_table', '--', 'LineWidth', 1);

% figure('Color', 'w')
% plot(time_table, cc_table, 'LineWidth', 3); hold on;




