function est = iestlsq(model, x, y)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Return structure which contains a function which estimates the input data
%   by least squares method.
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   x, y ... input data
%   model ... model of function which estimates the data
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   est.f ... function (can be interval function)
%   est.parameters ... list of parameters of function (estimator)
%   est.f_upper ... real function, upper boundry of the estimator 
%   est.f_lower ... real function, lower boundry of the estimator
%   est.f_string ... function printed in a string
%   est.f_upper_string ... upper function printed in a string
%   est.f_lower_string ... lower function printed in a string
%
%------------------------------------------------------------------------
% .Implementation details.
%
%   Used formulation of A. Neumaier, Linear interval equations (1986).
%
%   |Im A|.|y|  =  |b|
%   |A' 0| |x|     |0|
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
%   2016-10-16    first version
%   2016-11-16    function accept interval input
%   2017-08-06    take linearization part out of function
%
%------------------------------------------------------------------------
% .Todo.
%
%ENDDOC===================================================================

if (nargin != 3)
  print_usage();
  return;
endif

[A,b] = ilin(model, x, y);

[m, n] = size(A);

C = [  eye(m) A; A' zeros(n) ];
d = [ b; zeros( n,1) ];

x = C \ d;
parameters = x((m+1:m+n));

est = idelin(model, parameters);

endfunction
