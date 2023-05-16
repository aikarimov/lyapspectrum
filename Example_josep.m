%Example of lyapspectrum for the Josephson system
close all
T = 400;
h = 0.01;
tspan = 0:h:T-h;

rng(4);

y0 = [(6*rand - 3), (8*rand - 2), (2*rand)]';

[t,y] = ode78(@josep,tspan,y0);

y(:,1) = sin(y(:,1));

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
%y0 = y(end,:); %set start point as the last point
divisF = 10; %division factor
h = h*divisF;
tspan = 0:h:T-h;

[L,~,Lexp] = lyapspectrum(@josep,tspan,y0,'disp','3d','view',[1 0 0],'jacobian',@Jjosep,'df',divisF);

for i = 1:3
    disp([num2str(L(i)),'+/-',num2str(2*std(Lexp(i,:)))]);
end

figure;
plot(tspan,Lexp');
grid;
set(gca,'TickLabelInterpreter','latex');
xlabel('$t$, sec','interpreter','latex'); 
ylabel('$\lambda$','interpreter','latex');
legend('$\lambda_1$','$\lambda_2$','$\lambda_3$','interpreter','latex');