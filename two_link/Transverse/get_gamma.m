function gamma = get_gamma(in1)
%GET_GAMMA
%    GAMMA = GET_GAMMA(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    22-Apr-2021 14:52:49

phi1 = in1(1,:);
phi2 = in1(2,:);
gamma = sin(phi1+phi2).*(4.9e1./2.0);
