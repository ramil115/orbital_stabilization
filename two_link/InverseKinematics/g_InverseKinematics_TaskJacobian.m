function TaskJacobian = g_InverseKinematics_TaskJacobian(in1)
%G_INVERSEKINEMATICS_TASKJACOBIAN
%    TASKJACOBIAN = G_INVERSEKINEMATICS_TASKJACOBIAN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    31-Mar-2021 12:05:44

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
t2 = q1+q2+q3;
t3 = sin(t2);
t4 = q1+q2;
t5 = sin(t4);
TaskJacobian = reshape([1.0,0.0,t3.*(-1.0./2.0)-t5./2.0-sin(q1)./2.0,0.0,1.0,t3.*(-1.0./2.0)-t5./2.0,0.0,0.0,t3.*(-1.0./2.0)],[3,3]);