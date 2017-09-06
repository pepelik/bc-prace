function [x, y] = ilin(model, x, y)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Take input vectors and return linearized data which can be estimated by any
%   linear interval estimation. See function idelin() for delinearization of
%   parameters of estimator (it is the invers proces of the linearization).
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   x, y ... input data
%   model ... model which is used for transformation of the input data
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   x, y ... transformed data which can be estimated by linear estimation
%
%------------------------------------------------------------------------
% .Implementation details.
%
%   Implemented models are from Catalgoue of curves for curve fitting by Vera
%   Sit and Melanie Poulin-Costello.
%
%------------------------------------------------------------------------
%  .Licence.
%
%   Copyright (C) 2015  Charles University in Prague, Czech Republic
%
%   LIME 1.0 is free for private use and for purely academic purposes.
%   It would be very kind from the future user of LIME 1.0 to give
%   reference that this software package has been developed
%   at Charles University, Czech Republic.
%
%   For any other use of LIME 1.0 a license is required.
%
%   THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, WITHOUT LIMITATIONS, THE IMPLIED WARRANTIES
%   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%
%-------------------------------------------------------------------------
% .History.
%
%   2017-08-05    first version
%
%------------------------------------------------------------------------
% .Todo.
%
%ENDDOC===================================================================

if (nargin != 3)
  print_usage();
  return;
endif


[m,n] = size(x);

switch (model)
  case {'poly1', 'lin'} % y = p_1 + p_2*x
    x = [ones(m,1), x];
  case {'poly2', 'quad'} % y = p_1 + p_2*x + p_3*x^2
    x = [ones(m,1), x, x.^2];
  case {'poly3', 'cub'} % y = p_1 + p_2*x + p_3*x^2 + p_4*x^3
    x = [ones(m,1), x, x.^2, x.^3];
  case {'exptype1', 'exp'} %y = p_1*exp(p_2*x)
    x = [ones(m,1), x];
    y = log(y);
  case 'exptype2' % y = exp(p_1 - p_2*x)
    x = [ones(m,1), x];
    y = log(y);
  case 'exptype3' % y = p_1*exp(p_2/x)
    x = [ones(m,1), 1./x];
    y = log(y);
  case 'exptype4' % y = p_1*p_2^x
    x = [ones(m,1), x];
    y = log(y);
  case 'pow' % y = p_1*x^p_2
    x = [ones(m,1), log(x)];
    y = log(y);
  case 'exppowtype1' % y = p_1 * x^p_2 * p_3^x
    x = [ones(m,1), log(x), x];
    y = log(y);
  case 'exppowtype2' % y = p_1 * x^p_2 * exp(p_3 * x)
    x = [ones(m,1), log(x), x];
    y = log(y);
  case 'log' % y = p_1 + p_2 * ln(x)
    x = [ones(m,1), log(x)];
  case 'loglin' % y = p_1 + p_2*x + p_3*ln(x)
    x = [ones(m,1), x, log(x)];
  case 'explin' % y = p_1 + p_2*x + p_3*exp(x)
    x = [ones(m,1), x, exp(x)];
  case 'explintype2' % y = p_1 + p_2*x + p_3*exp(-x)
    x = [ones(m,1), x, exp(-x)];
  otherwise 
    error ("estimationlsq: unknown model of function");
endswitch

endfunction
