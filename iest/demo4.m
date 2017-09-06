function demo4()
%%%% Demonstration of outer estimation %%%%

disp("Demo4 started");
disp("Demonstration of outer estimation")
disp("Model 'loglin' used")

% load data
load dataset2.txt

% used model
model = "loglin";

% cut the end of data
x = x(1:32);
y = y(1:32);

%% plot the interval input data (x is real, y is interval)
%figure()
%plot(x,y)
%title("Dataset2");
%xlabel("Number of breath");
%ylabel("N_2");



%%% LSQ %%%
disp("\n=== LSQ ===")
est = iestlsq(model, x, y);
disp(est.f_string);
figure()
plot(x,y)
title("LSQ");
xlabel("Number of breath");
ylabel("N_2");
hold on;
fplot(est.f_upper, axis(),'r')
fplot(est.f_lower, axis(), 'r')
axis("auto"); % autoscale axis
%print -dsvg apl-nlin-lsq.svg

%%% Tolerance approach %%%
disp("\n=== Tolerance approach ===")
esttol = iesttol(model,x,y); 
disp(esttol.f_string)

figure()
hold on;
plot(x,y)
title("Terance approach");
xlabel("Number of breath");
ylabel("N_2");
fplot(esttol.f_upper, axis(),'r')
fplot(esttol.f_lower, axis(), 'r')
axis("auto"); % autoscale axis

%%%% Outer estimation (symetrical estimation)%%%
disp("\n=== Outer estimation  (symetrical estimation)===")
est = iestout(model, x, y, [1;0]);
disp(est.f_string);

figure()
plot(x,y)
title("Outer estimation (symetrical estimation)");
xlabel("Number of breath");
ylabel("N_2");
hold on;
fplot(est.f_upper, axis(),'r')
fplot(est.f_lower, axis(), 'r')
axis("auto"); % autoscale axis
%print -dsvg apl-nlin-lsq.svg

%%%% Outer estimation (bounds separately)%%%
disp("\n=== Outer estimation (bounds separately) ===")
est = iestoutsep(model, x, y, [1;0]);
disp(est.f_string);

figure()
plot(x,y)
title("Outer estimation (boudns separately)");
xlabel("Number of breath");
ylabel("N_2");
hold on;
fplot(est.f_upper, axis(),'r')
fplot(est.f_lower, axis(), 'r')
axis("auto"); % autoscale axis
%print -dsvg apl-nlin-lsq.svg

endfunction 
