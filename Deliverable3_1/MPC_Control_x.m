classdef MPC_Control_x < MPC_Control
  
  methods
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets (Ignore this before Todo 3.2)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % SET THE HORIZON HERE
      N = 30; % Determined empirically
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      con = [];
      obj = 0;
      
      M = [1; -1]; 
      m = [0.3; 0.3];
      
      F = [0 1 0 0; 0 -1 0 0]; 
      f = [0.035; 0.035];
      
      A = mpc.A; 
      B = mpc.B;

      Q = diag([0.5 5 5 20]);
      R = 0.01;
      
      [K, P, ~] = dlqr(A,B,Q,R);
      
      % MATLAB defines K as -K, so invert its signal
      K = -K;
      
      Xf = polytope([F; M*K],[f; m]);
      Acl = A+B*K;
      while 1
        prevXf = Xf;
        [T,t] = double(Xf);
        preXf = polytope(T*Acl,t);
        Xf = intersect(Xf, preXf);
        if isequal(prevXf, Xf)
            break
        end
      end
      
      [Ff,ff] = double(Xf);
      
      % Plotting the terminal sets
      figure;
      plot(Xf.projection(1:2), 'b');
      title('X_f for subsystem x');
      xlabel('$\dot{\beta}$', 'interpreter', 'latex');
      ylabel('\beta');
      
      figure;
      plot(Xf.projection(2:3), 'b');
      title('X_f for subsystem x');
      xlabel('\beta');
      ylabel('$\dot{x}$', 'interpreter', 'latex');
      
      figure;
      plot(Xf.projection(3:4), 'b');
      title('X_f for subsystem x');
      xlabel('$\dot{x}$', 'interpreter', 'latex');
      ylabel('x');
      
      for i=1:N-1
          con = con + (x(:,i+1) == A*x(:,i) + B*u(:,i));
          con = con + (F*x(:,i) <= f) + (M*u(:,i) <= m);
          obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
      end
      
      con = con + (Ff*x(:,N) <= ff);
      obj = obj + x(:,N)'*P*x(:,N);
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position (Ignore this before Todo 3.2)
      ref = sdpvar;            
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      con = [];
      obj = 0;
      
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'),...
          ref, {xs, us});
      
    end
  end
end
