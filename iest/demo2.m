function demo2()
%%%%% Demonstration of models %%%
%Try to find model for medical data
%Plot only LSQ model from real data set

disp("Demo2 started");
disp("Real estimation - fitting middles of dataset2.txt by LSQ\n")

% load data
load dataset2.txt
y = mid(y);


%grid for plotting
r = 3;
c = 3;
%all functions
models = ["lin"; "quad"; "cub"; "exp"; "exptype2"; "exptype3"; "exptype4"; "pow"; "exppowtype1"; "exppowtype2"; "log"; "loglin";  "explintype2" ];
%relevant functions
%models = ["quad"; "exp"; "exptype2"; "exptype3"; "exptype4"; "pow"; "exppowtype1"; "exppowtype2"; "log"; "loglin"; "explintype2"];
hax = [];
figure(1)
for i = 1 : r*c
	if (i < length(models))
		hax(i) = subplot(r,c,i);
		hold on;
		model = deblank(models(i,:));
		title(sprintf("Model '%s'", model));
		%scatter(x, y, 1, 'b', "filled");
		plot(x,y, '.');
		est = iestlsq(model, x,y);
		fplot(est.f,axis(), 'r');
		disp(model)
		disp(est.f_string);
		legend("off");
		%axis("auto"); % autoscale axis
	end
end
%print -dsvg demo2-1.svg %save picture

figure(2);
hold on;
for i = 1 : r*c
	if (i+(r*c) <= length(models))
		hax(i) = subplot(r,c,i);
		hold on;
		model = deblank(models(i+(r*c),:));
		title(sprintf("Model '%s'", model));
		%scatter(x, y, 1, 'b', "filled");
		plot(x,y, '.');
		est = iestlsq(model, x,y);
		fplot(est.f,axis(), 'r');
		disp(model)
		disp(est.f_string);
		legend("off");
		%axis("auto"); % autoscale axis
	end
end
%print -dsvg demo2-2.svg %save picture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters = estimationrlsq(model,x,y) % accept only real values and return real parameters
	[m,n] = size(x);

printf("%s\t\t",model);
	
switch (model)
  case{'poly1','lin'}
	X = [ones(m,1), x];
	parameters = X \ y;
	p = numtostr(parameters);
	printf("y = %s + %s*x\n", char(p(1)), char(p(2)));
  case {'poly2', 'quad'} % y = p_1 + p_2*x + p_3*x^2
    X = [ones(m,1), x, x.^2];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s*x + %s*x^2\n", char(p(1)), char(p(2)), char(p(3)));
  case {'poly3', 'cub'} % y = p_1 + p_2*x + p_3*x^2 + p_4*x^3
    X = [ones(m,1), x, x.^2, x.^3];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s*x + %s*x^2 + %s*x^3\n", char(p(1)), char(p(2)), char(p(3)), char(p(4)));
  case {'exptype1', 'exp'} %y = p_1*exp(p_2*x)
    X = [ones(m,1), x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [exp(parameters_tmp(1)); parameters_tmp(2)];
    p = numtostr(parameters);
    printf("y = %s*exp(%s*x)\n", char(p(1)), char(p(2)));
  case 'exptype2' % y = exp(p_1 - p_2*x)
    X = [ones(m,1), x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [parameters_tmp(1); -parameters_tmp(2)];
    p = numtostr(parameters);
    printf("y = exp(%s - %s*x)\n", char(p(1)), char(p(2)));
  case 'exptype3' % y = p_1*exp(p_2/x)
    X = [ones(m,1), 1./x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [exp(parameters_tmp(1)); parameters_tmp(2)];
    p = numtostr(parameters);
    printf("y = %s*exp(%s/x)\n", char(p(1)), char(p(2)));
  case 'exptype4' % y = p_1*p_2^x
    X = [ones(m,1), x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = exp(parameters_tmp);
    p = numtostr(parameters);
    printf("y = %s*%s^x\n", char(p(1)), char(p(2)));
  case 'pow' % y = p_1*x^p_2
    X = [ones(m,1), log(x)];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [exp(parameters_tmp(1)); parameters_tmp(2)];
    p = numtostr(parameters);
    printf("y = %s * x^%s\n", char(p(1)), char(p(2)));
  case 'exppowtype1' % y = p_1 * x^p_2 * p_3^x
    X = [ones(m,1), log(x), x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [exp(parameters_tmp(1)); parameters_tmp(2); exp(parameters_tmp(3))];
    p = numtostr(parameters);
    printf("y = %s * x^%s * %s^x \n", char(p(1)), char(p(2)), char(p(3)));
  case 'exppowtype2' % y = p_1 * x^p_2 * exp(p_3 * x)
    X = [ones(m,1), log(x), x];
    y = log(y);
    parameters_tmp = X \ y;
    parameters = [exp(parameters_tmp(1)); parameters_tmp(2); parameters_tmp(3)];
    p = numtostr(parameters);
    printf("y = %s * x^%s * exp(%s*x)\n", char(p(1)), char(p(2)), char(p(3)));
  case 'log' % y = p_1 + p_2 * ln(x)
    X = [ones(m,1), log(x)];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s * log(x)\n", char(p(1)), char(p(2)));
  case 'loglin' % y = p_1 + p_2*x + p_3*ln(x)
    X = [ones(m,1), x, log(x)];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s*x + %slog(x)\n", char(p(1)), char(p(2)), char(p(3)));
  case 'explin' % y = p_1 + p_2*x + p_3*exp(x)
    X = [ones(m,1), x, exp(x)];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s*x + %s*exp(x)\n", char(p(1)), char(p(2)), char(p(3)));
  case 'explintype2' % y = p_1 + p_2*x + p_3*exp(-x)
    X = [ones(m,1), x, exp(-x)];
    parameters = X \ y;
    p = numtostr(parameters);
    printf("y = %s + %s*x + %s*exp(-x)\n", char(p(1)), char(p(2)), char(p(3)));
  otherwise 
    error ("estimationlsq: unknown model of function");
endswitch

endfunction 


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

endfunction %end of function find_model
