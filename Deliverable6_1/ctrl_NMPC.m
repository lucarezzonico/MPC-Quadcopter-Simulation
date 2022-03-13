function [ctrl, traj] = ctrl_NMPC(quad)

    import casadi.*
    opti = casadi.Opti();
    
    obj=0;

    T_s = 1/5;
    N = 20;
    
    % −−−− decision variables −−−−−−−−−
    X = opti.variable(12,N+1); % state trajectory variables
    U = opti.variable(4, N); % control trajectory (throttle, brake)
    
    X0 = opti.parameter(12,1); % initial state
    REF = opti.parameter(4,1); % reference position [x,y,z,yaw]
    
    Q = diag([1 1 1 1 1 5 1 1 1 5 5 5000]);
    R = diag([0.0001 0.0001 0.0001 0.0001]);
    
    M = [0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 1;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         1 0 0 0;
         0 1 0 0;
         0 0 1 0];
    
    f_discrete = @(x,u) RK4(x,u,T_s,quad);
    
    opti.subject_to(X(:, 1) == X0);
    opti.subject_to(0 <= U <= 1.5);
    
    for i=1:N
      opti.subject_to(X(:, i+1) == f_discrete(X(:, i),U(:, i)));
      obj = obj + (X(:, i) - M * REF)'*Q*(X(:, i) - M * REF) + U(:, i)'*R*U(:, i);
    end
    
    opti.minimize(obj + (X(:, N+1) - M * REF)'*Q*(X(:,N+1) - M * REF));
    
    ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end

function u = eval_ctrl(x, ref, opti, X0, REF, X, U)

    % −−−− Set the initial state and reference −−−−
    opti.set_value(X0, x);
    opti.set_value(REF, ref);
    
    % −−−− Setup solver NLP −−−−−−
    ops = struct('ipopt', struct('print_level',0, 'tol', 1e-3), 'print_time', false);
    opti.solver('ipopt', ops);
    
    % −−−− Solve the optimization problem −−−− 
    sol = opti.solve();
    assert(sol.stats.success == 1, 'Error computing optimal input');
    
    u = opti.value(U(:,1));
    
    % Use the current solution to speed up the next optimization
    opti.set_initial(sol.value_variables());
    opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end

%Taken from 7.1 exercise
function [x_next] = RK4(X,U,h,quad)
   k1 = quad.f(X, U);
   k2 = quad.f(X+h/2*k1, U);
   k3 = quad.f(X+h/2*k2, U);
   k4 = quad.f(X+h*k3, U);
   x_next = X + h/6*(k1+2*k2+2*k3+k4);
end