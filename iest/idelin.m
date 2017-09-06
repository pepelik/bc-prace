function retval = idelin(model, parameters)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Return structure which contains list of parameters of function which
%   estimates the input data. It is inverse operation of linearization by
%   function ilin().
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   model ... model used for linearization of the function which estimates the data
%   parameters ... parameters of linearized function
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   retval.param ... list of parameters of estimator (the delinearized function)
%   retval.f ... function (can be interval function)
%   retval.f_upper ... real function, upper boundry of the retval.f
%   retval.f_lower ... real function, lower boundry of the retval.f
%   retval.f_string ... function printed in a string
%   retval.f_upper_string ... upper function printed in a string
%   retval.f_lower_string ... lower function printed in a string
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

if (nargin != 2)
  print_usage();
  return;
endif

if (isa(parameters,"infsup"))

	fce = retf(model, parameters);
	fce_lower = retf(model, inf(parameters));
	fce_upper = retf(model, sup(parameters));
else
	fce = retf(model, parameters);
	fce_lower = fce;
	fce_upper = fce;
end

	retval = struct('param', fce.param, 'f', fce.f, 'f_string', fce.f_string, 
	                'f_upper', fce_upper.f, 'f_upper_string', fce_upper.f_string, 'f_upper_param', fce_upper.param, 
					'f_lower', fce_lower.f, 'f_lower_string', fce_lower.f_string, 'f_lower_param', fce_lower.param);

endfunction

function retval = retf(model, parameters)
	switch (model)
	  case {'poly1', 'lin'} % y = p_1 + p_2*x
		param = parameters;
		f = @(x) param(1) + param(2)*x;
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x\n", char(p(1)), char(p(2)));
	  case {'poly2', 'quad'} % y = p_1 + p_2*x + p_3*x^2
		param = parameters;
		f = @(x) param(1) + param(2)*x + param(3)*x.^2;
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x + %s*x^2\n", char(p(1)), char(p(2)), char(p(3)));
	  case {'poly3', 'cub'} % y = p_1 + p_2*x + p_3*x^2 + p_4*x^3
		param = parameters;
		f = @(x) param(1) + param(2)*x + param(3)*x.^2 + param(4)*x.^3;
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x + %s*x^2 + %s*x^3\n", char(p(1)), char(p(2)), char(p(3)), char(p(4)));
	  case {'exptype1', 'exp'} %y = p_1*exp(p_2*x)
		param = [exp(parameters(1)); parameters(2)];
		f = @(x) param(1)*exp(param(2)*x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s*exp(%s*x)\n", char(p(1)), char(p(2)));
	  case 'exptype2' % y = exp(p_1 - p_2*x)
		param = [parameters(1); -parameters(2)];
		f = @(x) exp(param(1) - param(2).*x);
		p = numtostr(parameters);
		f_string = sprintf("y = exp(%s - %s*x)\n", char(p(1)), char(p(2)));
	  case 'exptype3' % y = p_1*exp(p_2/x)
		param = [exp(parameters(1)); parameters(2)];
		f = @(x)  param(1)*exp(param(2)./x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s*exp(%s/x)\n", char(p(1)), char(p(2)));
	  case 'exptype4' % y = p_1*p_2^x
		param = exp(parameters);
		f = @(x) param(1)*param(2).^x;
		p = numtostr(parameters);
		f_string = sprintf("y = %s*%s^x\n", char(p(1)), char(p(2)));
	  case 'pow' % y = p_1*x^p_2
		param = [exp(parameters(1)); parameters(2)];
		f = @(x) param(1)*x.^param(2);
		p = numtostr(parameters);
		f_string = sprintf("y = %s * x^%s\n", char(p(1)), char(p(2)));
	  case 'exppowtype1' % y = p_1 * x^p_2 * p_3^x
		param = [exp(parameters(1)); parameters(2); exp(parameters(3))];
		fexppowtype1= @(x) param(1) * x^param(2) * param(3)^x;
		f = @(x) arrayfun(fexppowtype1,x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s * x^%s * %s^x \n", char(p(1)), char(p(2)), char(p(3)));
	  case 'exppowtype2' % y = p_1 * x^p_2 * exp(p_3 * x)
		param = [exp(parameters(1)); parameters(2); parameters(3)];
		fexppowtype2 = @(x) param(1) * x^param(2) * exp(param(3) * x);
		f = @(x) arrayfun(fexppowtype2, x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s * x^%s * exp(%s*x)\n", char(p(1)), char(p(2)), char(p(3)));
	  case 'log' % y = p_1 + p_2 * ln(x)
		param = parameters;
		f = @(x)  param(1) + param(2) * log(x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s * log(x)\n", char(p(1)), char(p(2)));
	  case 'loglin' % y = p_1 + p_2*x + p_3*ln(x)
		param = parameters;
		f = @(x) param(1) + param(2)*x + param(3)*log(x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x + %slog(x)\n", char(p(1)), char(p(2)), char(p(3)));
	  case 'explin' % y = p_1 + p_2*x + p_3*exp(x)
		param = parameters;
		f = @(x) param(1) + param(2)*x + param(3)*exp(x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x + %s*exp(x)\n", char(p(1)), char(p(2)), char(p(3)));
	  case 'explintype2' % y = p_1 + p_2*x + p_3*exp(-x)
		param = parameters;
		f = @(x) param(1) + param(2)*x + param(3)*exp(-x);
		p = numtostr(parameters);
		f_string = sprintf("y = %s + %s*x + %s*exp(-x)\n", char(p(1)), char(p(2)), char(p(3)));
	  otherwise 
		error ("estimationlsq: unknown model of function");
	endswitch

	retval = struct('param', param, 'f', f, 'f_string', f_string);
endfunction %retf


function s = numtostr(v)
% convert number vector (real or interval) to cellstr
	s = cellstr([]);
	if (isa(v, "infsup"))
	  for i = 1:length(v)
		s(i) = sprintf("[%d, %d]", inf(v(i)), sup(v(i)));
	  end
	  else
	  for i = 1:length(v)
		s(i) = num2str(v(i));
	  end
	endif
endfunction
