function est = iestoutsep(model, x, y, tol)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Returned structure includes functions which are bounds of interval curve.
%   The curve includes all input data (and the bounds are founded separately).
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   model ... model of estimation
%   x, y  ... vectors of input data (must be the same size) 
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
%   est.est_real ... the central real estimator
%   est.y_lin ... shifted input data
%
%------------------------------------------------------------------------
% .Implementation details.
%
%   The function estimate centers of data by least square method.
%   Then linearize the data (diference between y and estimator) and estimate it
%   by tolerance approach (function linregtol) - but every boundary function
%   independently.  Returned function is defined by upper and lower bounds. It
%   is sum of the estimators.
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
%   2017-06-15   first version 
%   2017-08-09   return structure
%
%------------------------------------------------------------------------
% .Todo.
%   Interval data case.
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

% two subproblems (one for upper boundary, second for lower boundary), both are real problems
y_up = [];
x_up = [];
y_low = [];
x_low = [];
for i = flip([1:length(x)]) % start at the end and split the data by the real estimator
	if (isa(y,"infsup"))
		if (sup(y_lin(i)) > 0) % upper bound of interval is positive
			y_up = [sup(y_lin(i)); y_up];
			x_up = [x(i); x_up];
		endif
		if (inf(y_lin(i)) < 0) % lower bound is negative
			y_low = [inf(y_lin(i)); y_low];
			x_low = [x(i); x_low];
		endif
	  else %y is real vector
		if (y_lin(i) > 0)
			y_up = [y_lin(i); y_up];
			x_up = [x(i); x_up];
		endif
		if (y_lin(i) < 0)
			y_low = [y_lin(i); y_low];
			x_low = [x(i); x_low];
		endif
	endif
endfor

if (isa(x, "infsup"))
	y_up = [y_up ; y_up];
	x_up = [inf(x_up); sup(x_up)];

	y_low = [y_low ; y_low];
	x_low = [inf(x_low); sup(x_low)];
endif

% solving the two linear real problems
if (nargin < 4)
	lin_est_up = iesttol('lin', x_up, y_up);
	lin_est_low = iesttol('lin', x_low, y_low);
else 
	lin_est_up = iesttol('lin', x_up, y_up,tol);
	lin_est_low = iesttol('lin', x_low, y_low,tol);
endif

% delinearization
f_upper = @(x) est.f(x) + lin_est_up.f_upper(x);
f_lower = @(x) est.f(x) + lin_est_low.f_lower(x);

% construction of output structure (sum of central crisp estimation function and linear estimator)
f = @(x) infsup(f_lower(x), f_upper(x));
f_string = sprintf("\ny_upper = g(x) + h(x)\n g(x) = %s\n h(x) = %s\ny_lower = g(x) + h(x)\n g(x) = %s\n h(x) = %s", 
           est.f_string, lin_est_up.f_upper_string, 
		   est.f_string, lin_est_low.f_lower_string);

f_upper_string = sprintf("y_upper = g(x) + h(x)\n g(x) = %s\n h(x) = %s", est.f_string, lin_est_up.f_upper_string);
f_lower_string = sprintf("y_lower = g(x) + h(x)\n g(x) = %s\n h(x) = %s", est.f_string, lin_est_low.f_lower_string);

est = struct('f', f, 'f_string', f_string,
             'f_upper', f_upper, 'f_upper_string', f_upper_string,
			 'f_lower', f_lower, 'f_lower_string', f_lower_string,
			 'est_real', est,
			 'y_lin', y_lin);

endfunction
