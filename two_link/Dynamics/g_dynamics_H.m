function H = g_dynamics_H(in1)
%G_DYNAMICS_H
%    H = G_DYNAMICS_H(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    22-Apr-2021 14:52:48

q2 = in1(2,:);
t2 = cos(q2);
t3 = t2.*(5.0./4.0);
t4 = t3+2.5e1./2.4e1;
H = reshape([t2.*(5.0./2.0)+5.5e1./1.2e1,t4,t4,2.5e1./2.4e1],[2,2]);
