function est = iestout(model, x, y, tol)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Returned structure includes interval function which bounds define interval
%   curve. The curve includes all input data. It is the outer estimator of the 
%   dataset.
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   model ... model of estimation
%   x, y  ... vectors of input data 
%	tol   ... tolerance rate vector for linregtol (relative tolerance is 
%             used if it is not given)
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   est.f ... interval function which is interval outer estimator of the data
%   est.f_upper ... upper bound
%   est.f_lower ... lower bound
%   est.f_string ... function printed in a string
%   est.f_upper_string ... upper function printed in a string
%   est.f_lower_string ... lower function printed in a string
%   est.param ... parameters of function
%   est.est_real ... the central real estimator
%   est.est_lin ... the linear interval estimator of shifted data
%   est.y_lin ... shifted input data
%
%------------------------------------------------------------------------
% .Implementation details.
%
%   The function estimate centers of data by least square method.
%   Then linearize the data (difference between y and estimator) and estimate it
%   by tolerance approach (function linregtol). Returned structure includes
%   function it is a sum of the estimators.
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
%   2016-11-17    first version
%   2017-06-15    renamed function (from nlinregtol) and added argument tol
%   2017-08-09    return structure not only interval function 
%
%------------------------------------------------------------------------
% .Todo.
%
%ENDDOC===================================================================

if (nargin > 4)
	print_usage();
	return;
end

% lineariation
if (isa(x,"infsup"))
	mid_x = mid(x);
  else
  	mid_x = x;
end

if (isa(y,"infsup"))
	mid_y = mid(y);
  else
  	mid_y = y;
end

est = iestlsq(model, mid_x, mid_y);

y_lin = y - est.f(x);

% solution of the linear problem
if (nargin < 4)
	lin_est = iesttol('lin',x, y_lin);
else 
	lin_est = iesttol('lin',x, y_lin, tol);
endif


% sum of central crisp estimation function and linear interval function
f = @(x) est.f(x) + lin_est.f(x);
f_upper = @(x) est.f(x) + lin_est.f_upper(x);
f_lower = @(x) est.f(x) + lin_est.f_lower(x);

f_string = sprintf("y = g(x) + h(x)\n g(x) = %s\n h(x) = %s", est.f_string, lin_est.f_string);
f_upper_string = sprintf("y = g(x) + h(x)\n g(x) = %s\n h(x) = %s", est.f_string, lin_est.f_lower_string);
f_lower_string = sprintf("y = g(x) + h(x)\n g(x) = %s\n h(x) = %s", est.f_string, lin_est.f_lower_string);

param = [est.param ; lin_est.param];
est_real = est;
est_lin = lin_est;

est = struct('f', f, 'f_upper', f_upper, 'f_lower', f_lower, 'param', param,
			 'f_string', f_string, 'f_upper_string', f_upper_string, 'f_lower_string', f_lower_string,
			 'est_real', est_real,
			 'est_lin', est_lin,
			 'y_lin', y_lin);

endfunction
