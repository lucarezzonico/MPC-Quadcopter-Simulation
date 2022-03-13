classdef MPC_Control_yaw < MPC_Control
  
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
      N = 30;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      con = [];
      obj = 0;
      
      M = [1; -1]; 
      m = [0.2; 0.2];
      
      A = mpc.A; 
      B = mpc.B;

      Q = diag([0.5 5]);
      R = 0.1;
      
      [K,P,~] = dlqr(A,B,Q,R);
      
        
      % MATLAB defines K as -K, so invert its signal
      K = -K;
      
      Xf = polytope(M*K, m);
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
      
      for i=1:N-1
          con = con + ((x(:,i+1) - xs) == A*(x(:,i) - xs) + B*(u(:,i) - us));
          con = con + (M*u(:,i) <= m);
          obj = obj + (x(:,i) - xs)'*Q*(x(:,i) - xs) + (u(:,i) - us)'*R*(u(:,i) - us);
      end
      
      con = con + (Ff*x(:,N) <= ff + Ff*xs);
      obj = obj + (x(:,N) - xs)'*P*(x(:,N) - xs);
      
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
            
      M = [1;-1]; 
      m = [0.2;0.2];
      
      A = mpc.A; 
      B = mpc.B;
      C = mpc.C; 
      D = mpc.D;
      
      con = [];
      obj = 0;
      
      con = con + (xs == A*xs + B*us);
      con = con + (ref == C*xs + D*us);
      con = con + (M*us <= m);
      
      obj = obj + us^2;
      
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'),...
          ref, {xs, us});
      
    end
  end
end
