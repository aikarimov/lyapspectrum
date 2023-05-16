%Example of lyapspectrum for the GRN system
close all
T = 200;
h = 0.01;
tspan = 0:h:T-h;

y0 = [0.4, 0.39, 0.4, 0.395]';

[t,y] = ode78(@GRN,tspan,y0);

figure(1);
plot(t,y);
xlabel('$t$, sec','interpreter','latex'); ylabel('$x$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

figure(3);
plot3(y(:,1),y(:,2),y(:,3));
xlabel('$\sin(x_1)$','interpreter','latex'); 
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
view([1 0 0]);

%Lyapunov calculation
divisF = 10; %division factor
h = h*divisF;
tspan = 0:h:T-h;

[L,~,Lexp] = lyapspectrum(@GRN,tspan,y0,'disp','3d','view',[1 0 0],'df',divisF);

for i = 1:4
    disp([num2str(L(i)),'+/-',num2str(2*std(Lexp(i,:)))]);
end

figure(10);
semilogy(tspan,Lexp(1,:));

figure(11);
plot(tspan,Lexp');
xlabel('$t$, sec','interpreter','latex');
ylabel('$\lambda_i$','interpreter','latex');
grid;