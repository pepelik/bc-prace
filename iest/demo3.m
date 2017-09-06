%%%%% Demonstration of interval estimation LSQ%%%
%
%Plot interval LSQ models

disp("Demo3 started");
disp("Ordinary least squares estimation")

% load data
load dataset2.txt

%grid for plot graphs
r = 3;
c = 3;
% models
models = ["lin"; "quad"; "cub"; "exp"; "exptype2"; "exptype3"; "exptype4"; "pow"; "exppowtype1"; "exppowtype2"; "log"; "loglin"; "explintype2"];

hax = [];
figure(1)
for i = 1 : r*c
	if (i < length(models))
		hax(i) = subplot(r,c,i);
		hold on;
		model = deblank(models(i,:));
		title(sprintf("Model '%s'", model));
		plot(x,y);
		est = iestlsq(model, x,y);
		fplot(est.f_upper, axis(),'r')
		fplot(est.f_lower, axis(), 'r')
		legend("off");
	end
end

figure(2);
hold on;
for i = 1 : r*c
	if (i+(r*c) <= length(models))
		hax(i) = subplot(r,c,i);
		hold on;
		model = deblank(models(i+(r*c),:));
		title(sprintf("Model '%s'", model));
		plot(x,y);
		est = iestlsq(model, x,y);
		fplot(est.f_upper, axis(), 'r')
		fplot(est.f_lower, axis(), 'r')
		legend("off");
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
