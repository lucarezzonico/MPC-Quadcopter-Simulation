classdef MPC_Control_z < MPC_Control
  properties
    A_bar, B_bar, C_bar % Augmented system for disturbance rejection    
    L                   % Estimator gain for disturbance rejection
  end
  
  methods
    function mpc = MPC_Control_z(sys, Ts)
      mpc = mpc@MPC_Control(sys, Ts);
      
      [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
    end
    
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   d_est  - disturbance estimate
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets (Ignore this before Todo 3.3)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);

      % SET THE HORIZON HERE
      N = 30;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      
      con = [];
      obj = 0;
      
      M = [1; -1]; 
      m = [0.3; 0.2];
      
      A = mpc.A; 
      B = mpc.B;

      Q = diag([0.5 5]);
      R = 0.1;
      
      [~, P, ~] = dlqr(A,B,Q,R);
      
      for i=1:N-1
          con = con + ((x(:,i+1) - xs) == A*(x(:,i) - xs) + B*(u(:,i) - us));
          con = con + (M*u(:,i) <= m);
          obj = obj + (x(:,i) - xs)'*Q*(x(:,i) - xs) + (u(:,i) - us)'*R*(u(:,i) - us);
      end
 
      obj = obj + (x(:,N) - xs)'*P*(x(:,N) - xs); 
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us, d_est}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      %   d_est  - disturbance estimate
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position (Ignore this before Todo 3.3)
      ref = sdpvar;
            
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);
      
      con = [];
      obj = 0;
      
      A = mpc.A; 
      B = mpc.B;
      C = mpc.C; 
      D = mpc.D;
      
      M = [1; -1]; 
      m = [0.3; 0.2];
      
      con = con + (xs == A*xs + B*us);
      con = con + (ref == C*xs + D*us);
      con = con + (M*us <= m);
      
      obj = obj + us^2;
      
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'),...
          {ref, d_est}, {xs, us});
    end
    
    
    % Compute augmented system and estimator gain for input disturbance rejection
    function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
      
      %%% Design the matrices A_bar, B_bar, L, and C_bar
      %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
      %%% converges to the correct state and constant input disturbance
      %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      
      A_bar = [];
      B_bar = [];
      C_bar = [];
      L = [];
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
  end
end
