function h = get_handler_full_dnmcs_ode(HCg_model)
h = @(t,x,u) full_dnmcs_ode(x,u,HCg_model);
% h = @full_dnmcs_ode;
% function  dxdt = full_dnmcs_ode(x,u)
function  dxdt = full_dnmcs_ode(x,u,HCg_model)
    n = HCg_model.dof_configuration_space_robot;
    
    H = HCg_model.get_H(x(1:n));
    C = HCg_model.get_C(x(1:n), x((n+1):end));
    gr = HCg_model.get_g(x(1:n));        
    T = HCg_model.get_T();   
    dxdt = [x((n+1):end); H\(T * u - C * x((n+1):end) - gr)];
end
end