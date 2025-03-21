function dsys = cont2disc(sys, Ts, method)

% CONT2DISC
%
% dsys = cont2disc(sys, Ts, method)
%
% Discretize a continuous-time dynamic system. sys can be a transfer
% function or a state-space model and Ts is the sampling time. The
% third input arg method determines the discretization method:
%
% 'backward' - backward difference
% 'forward'  - forward difference
% 'euler'    - same as 'forward'
% 'tustin'   - tustin method
% 'bilinear' - same as 'tustin'
% 'impulse'  - impulse-invariance
% 'step'     - step-invariance
% 'ramp'     - ramp-invariance
% 'zoh'      - zero-order hold
%
% This function does roughly the same as Matlab's c2d but it adds a
% few methods. However, only SISO models with no delays work at the
% moment.
%
% Note that if you use the 'impulse' method, the results will be
% different from using Matlab's c2d because of how c2d scales the
% discretized transfer function (see 'help c2d').
%
% Needs the Control and Symbolic toolboxes.
  
% Copyright (C) 2023-2025 Toni Saarela
% 2023-10-25 - ts - written
% 2024-12-18 - ts - added impulse, step, and ramp invariance methods
% 2025-01-25 - ts - added zoh;
%                   switched order of inputs sys,method,Ts -> sys,Ts,method
% 2025-03-20 - ts - output ss model if input is ss model; wrote help
  
% TODO:
% - Add Tustin prewarping
% - Handle MIMO
% - Handle delays

  % Symbolic variables needed
  syms z s k t

  % Was the input a state-space model?
  isss = isa(sys, 'ss');
  
  % ZOH discretization is done using numerical integration, so skip
  % going symbolic in that case. Otherwise, convert the transfer
  % function to symblic form.
  if ~strcmp(method,'zoh')
    if isss
      sys = tf(sys);
    end
    ps = tf2sym(sys);
  end

  switch method
    case 'backward'
      % The backward/forward Euler and Tustin methods are done with
      % appropriate substitutions.
      zsub = (z-1)/(z*Ts);
      pz = subs(ps, s, zsub);
    case {'forward', 'euler'}
      zsub = (z-1)/(Ts);
      pz = subs(ps, s, zsub);
    case {'tustin','bilinear'}
      zsub = 2/Ts*(z-1)/(z+1);
      pz = subs(ps, s, zsub);
    case 'impulse'
      % In the impulse, step, and ramp invariance methods, go to time
      % domain with inverse Laplace, do the appropriate substitution,
      % then do z-transform (appropriately scaled for step and ramp).
      ht = ilaplace(ps);
      ksub = k*Ts;
      hk = subs(ht, t, ksub);
      % NOTE: Matlab's c2d function scales the discretized transfer
      % function by the sampling interval for some reason, maybe so
      % that Matlab's impulse function produces the correct impulse
      % response? To get that same behavior, you'd need to do:      
      % hk = subs(ht*Ts, t, ksub)
      pz = ztrans(hk);
    case 'step'
      ht = ilaplace(1/s*ps);
      ksub = k*Ts;
      hk = subs(ht, t, ksub);
      pz = (z-1) / z * ztrans(hk);
    case 'ramp'
      ht = ilaplace(ps / s^2);
      ksub = k*Ts;
      hk = subs(ht, t, ksub);
      pz = (z-1)^2 / z * ztrans(hk);
    case 'zoh'
      % Zero-order hold is handled differently, not with substitutions
      % as above
      sys = ss(sys);

      % The following computes the matrix exponential of A*Ts (Phi)
      % and the integral of that times B (Gamma)
      [ma,na] = size(sys.A);
      [mb,nb] = size(sys.B); % ma==mb...
      M = expm([[sys.A sys.B]*Ts; zeros([nb na+nb])]);
      Phi = M(1:na,1:na);
      Gamma = M(1:na,na+1:na+nb);
      
      dsys = ss(Phi, Gamma, sys.C, sys.D, Ts);

      % If input was a transfer function, return one 
      if ~isss
        dsys = tf(dsys);
      end
      
      % not using symbolic things here, so return right away
      return
    otherwise
      error('Unknown discretization method.');
  end

  % Convert the symbolic model to Matlab TF object
  dsys = sym2tf(pz, Ts);

  % If input was a state space model, return one
  if isss
    dsys = ss(dsys);
  end
  
end
