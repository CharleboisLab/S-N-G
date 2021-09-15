function [t, X] = ODE_solver_model1_scenario1(t_end,dt,S,N,G1)
%Solves a system of coupled ODEs via MATLAB solver (default solver: ode45)
[t X] = ode45(@equations,0:dt:t_end,[S N G1]); 

end


function dx = equations(t,x)
    
    global kS rSN rNS rG1S dS kN rG1N dN kG1 dG1 
    
    dx = zeros(3,1);
    
    dx(1) = x(1)*kS + x(2)*rSN - x(1)*rNS - x(1)*rG1S - x(1)*dS;
    dx(2) = x(2)*kN + x(1)*rNS - x(2)*rG1N - x(2)*rSN - x(2)*dN;
    dx(3) = x(3)*kG1 + x(1)*rG1S + x(2)*rG1N - x(3)*dG1; 
    
end

