%Example of lyapspectrum for the Lorenz system

T = 30;
h = 0.01;
tspan = 0:h:T-h;

y0 = [8 12 12]';

[L,~,Lexp] = lyapspectrum(@lorenz2,tspan,y0,'disp','all','jacobian',@Jlorenz2,'df',10);

for i = 1:3
    disp([num2str(L(i)),'+/-',num2str(2*std(Lexp(i,:)))]);
end