function sys = sym2tf(p, Ts)

% SYM2TF
%
% sys = sym2tf(p, Ts)
% 
% Convert a transfer function given in symbolic form to Matlab's
% transfer function. The variable should be s for continuous and z for
% discrete transfer functions (althoug any variable other than s will
% give a discrete tf).
  
% Copyright (C) 2023-2025 Toni Saarela
% 2023-10-25 - ts - written
% 2024-12-18 - ts - added minreal
% 2025-03-20 - ts - help

  iscont = symvar(p)=='s';

  % p = simplify(p);
  [snum, sden] = numden(p);

  num = sym2poly(snum);
  den = sym2poly(sden);

  if iscont
    sys = tf(num, den);
  else
    sys = tf(num, den, Ts);
  end

  sys = minreal(sys);
  
end
