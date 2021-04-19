function g = g_dynamics_g(in1)
%G_DYNAMICS_G
%    G = G_DYNAMICS_G(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    16-Apr-2021 12:01:35

q1 = in1(1,:);
q2 = in1(2,:);
t2 = q1+q2;
t3 = sin(t2);
t4 = t3.*(4.9e1./2.0);
g = [t4+sin(q1).*(1.47e2./2.0);t4];