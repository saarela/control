function [p, Ts] = tf2sym(sys)

% TF2SYM
%
% [p, Ts] = tf2sym(sys)
%
% Convert a transfer function to symbolic form. The variable will be s
% for continuous and z for discrete transfer functions.
    
% Copyright (C) 2023-2025 Toni Saarela
% 2023-10-25 - ts - written
% 2025-03-20 - ts - help
  
  iscont = sys.variable=='s';
  Ts = sys.Ts;
  
  % sys = minreal(sys);
  [num, den] = tfdata(sys, 'v');

  if iscont
    syms s
    snum = poly2sym(num, s);
    sden = poly2sym(den, s);
  else
    syms z
    snum = poly2sym(num, z);
    sden = poly2sym(den, z);
  end
  
  p = snum / sden;
  
end
