%%%%% Demonstration of linear estimation %%%
%
%Plot only LSQ model from real data set

disp("Demo1 started");
disp("Linear outer estimation")

% load data
load dataset1.txt;

%%% Tolerance approach %%%
disp("\n=== Tolerance approach ===")
plot(x,y, '.');
axis([290 315 1240 1380])
hold on;

est = iesttol('lin', x, y);
fplot(est.f_upper, axis(), 'r')
fplot(est.f_lower, axis(), 'r')
disp(est.f_string);
%print -dsvg obr-lin-est-tol.svg

%%%% Outer estimation (bounds separately)%%%
disp("\n=== Outer estimation, bounds separately ===")
estoutsep = iestoutsep('lin', x, y);
disp(estoutsep.f_string);

figure()
clf;
hold on;
plot(x,y,'.')
fplot(estoutsep.f_lower, axis(), 'r');
fplot(estoutsep.f_upper, axis(), 'r');
axis([290 315 1240 1380])
%print -dsvg obr-lin-est-outer.svg
