function est = iesttol(model, x, y, tol, est)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikánová
%
%-------------------------------------------------------------------------
% .Description.
%
%   Computes verified bounds on linear possibilistic (strong) regression model
%   X betas = y by tolerance approach.
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   x ... input data vector
%   y ... vector of output data
%   tol ... tolerance rate vector
%   est ... parameters of central estimator
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   est ... structure of return values
%   
%   est.f ... interval function
%   est.param ... list of parameters of estimator
%   est.f_upper ... real function, upper boundry of the estimator
%   est.f_lower ... real function, lower boundry of the estimator
%   est.f_string ... function printed in a string
%   est.f_upper_string ... upper function printed in a string
%   est.f_lower_string ... lower function printed in a string
%
%------------------------------------------------------------------------
% .Implementation details.
%
%  see Hladík, Milan, and Michal Černý. "Interval regression by tolerance
%  analysis approach." Fuzzy Sets and Systems 193 (2012): 85-107.
%
%   Possibilistic model is implemented.
%   For crisp input-interval ouptu is used (strong) possibilistic model.
%   For interval input-interval output is used weak possibilistict model.
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
%   2016-07-11    first version
%   2016-08-06    added models which can be linearized
%
%------------------------------------------------------------------------
% .Todo.
%
%   Error if interval input-crisp output are given
%
%ENDDOC===================================================================

[m, n] = size(x);
if (m ~= length(y) )
   error('x and y must be of the same size.');
end

% linearization
[X,y] = ilin(model,x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crisp input - crisp output
if (~isa(X,"infsup") && ~isa(y, "infsup"))
	% central line
	if (nargin < 5)
		est = iestlsq(model,x,y).param;
	end

	% tolerance vector (is absolute value of central vector if it is not given)
	if (nargin < 4)
		tol = abs(est);
	end
	
  	retval_tolcc = linregtolcc(X,y,est,tol);
	delta = retval_tolcc.maxdelta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crisp input - interval output
if (~isa(X,"infsup") && isa(y, "infsup"))
	% central line
	if (nargin < 5)
		est = iestlsq(model,x,mid(y)).param;
	end
	% tolerance vector (is absolute value of central vector if it is not given)
	if (nargin < 4)
		tol = abs(est);
	end
	
	for j = 1:m
		if abs(inf(y(j)) - X(j,:) * est) >= abs( sup(y(j)) - X(j, :) * est)
			y_crisp(j) = inf(y(j));
		else
			y_crisp(j) = sup(y(j));
		end
	end
	y_crisp = y_crisp';
	%crisp-crisp
  	retval_tolcc = linregtolcc(X,y_crisp,est,tol);
	delta = retval_tolcc.maxdelta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interval input - interval output
[m, n] = size(X);

if (isa(X,"infsup") && isa(y, "infsup"))
	if (nargin < 5)
		est = iestlsq(model, mid(x), mid(y)).param; % estimator (parameters of "the central line")
	end
	if (nargin < 4)
		tol = abs(est);
	end

	%1. reduction
	X1 = zeros(m,n);
	y1 = zeros(m);
	for j = 1:m
		y1 = inf(y);
		for i = 1:n
			if est(i) >= 0
				X1(j,i) = sup(X(j,i));
			else
				X1(j,i) = inf(X(j,i));
			end 
		end
	end


	%2. reduction
	X2 = zeros(m,n);
	y2 = zeros(m);
	for j = 1:m
		y2 = sup(y);
		for i = 1:n
			if est(i) >= 0
				X2(j,i) = inf(X(j,i));
			else
				X2(j,i) = sup(X(j,i));
			end 
		end
	end

  	retval_tolcc1 = linregtolcc(X1,y1,est,tol);
	delta1 = retval_tolcc1.maxdelta;
  	retval_tolcc2 = linregtolcc(X2,y2,est,tol);
	delta2 = retval_tolcc2.maxdelta;
	delta = max(delta1, delta2);
	
	deltas = max(retval_tolcc1.deltas, retval_tolcc2.deltas);
	retval_tolcc = struct('deltas', deltas);
end %end of i-i case

betas = infsup(est - delta*tol, est + delta*tol);

retval_idelin = idelin(model, betas);
est = struct('param', retval_idelin.param, 'f', retval_idelin.f, 'f_string', retval_idelin.f_string, 
         'f_upper', retval_idelin.f_upper, 'f_upper_string', retval_idelin.f_upper_string, 'f_upper_param', retval_idelin.f_upper_param,
		 'f_lower', retval_idelin.f_lower, 'f_lower_string', retval_idelin.f_lower_string, 'f_lower_param', retval_idelin.f_lower_param,
		 'deltas', retval_tolcc.deltas, 'maxdelta', delta);

endfunction % iesttol()


function retval = linregtolcc(X, y, c, tol)
% LINear REGression by TOLerance approach for Crisp input-Crisp output
%
% Computes delta constant on linear regression model
% X betas = y by tolerance approach, crisp input - crisp output.
%
% Input parameters
%   X ... crisp input data matrix
%   y ... crisp vector of output data
%   c ... central vector (define central regression line)
%   tol ... tolerance rate vector
%
% Output parameters
%   maxdelta ... parametr of radius, real number
%   deltas ... all deltas
%
%------------------------------------------------------------------------
	[m, n] = size(X);

	maxdelta = 0;
	deltas = [];
	for j = 1:m
	  if (abs(X(j, :))*tol > 0)
		  newdelta = ( abs(y(j) - ( X(j, :)*c )) ) / ( abs(X(j, :))*tol );
		  deltas(j) = newdelta;
		  if (newdelta > maxdelta)
			 maxdelta = newdelta;
		  end
	  end
	end
%	delta = maxdelta;
	retval = struct('maxdelta', maxdelta, 'deltas', deltas);
end


