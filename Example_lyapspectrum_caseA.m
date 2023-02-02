%Example of lyapspectrum for the Lorenz system

T = 300;
h = 0.1;
tspan = 0:h:T-h;

y0 = [0 5 0]';

[L,~,Lexp] = lyapspectrum(@CaseA,tspan,y0,'disp','all','view',[0 1 0],'jacobian',@JCaseA,'df',10);

for i = 1:3
    disp([num2str(L(i)),'+/-',num2str(2*std(Lexp(i,:)))]);
end