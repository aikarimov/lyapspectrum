function [L, Lspan, Lexp] = lyapspectrum(varargin)
%LYAPSPECTRUM Calculate the Lyapunov spectrum for a particular system
%   [L, LSPAN, LEXP]  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0)
%   [L, LSPAN, LEXP]  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0,'NAME',VALUE)
%   Input:
%   ODEFUN is a function handle ODEFUN(T,Y)
%   TSPAN = [T0 T1 ... TFINAL] is a vector of times
%   Y0 is a vector of initial conditions
%   'NAME',VALUE given an additional input options:
%   'jacobian',JFCN is a Jacobian matrix ODEFUN(T,Y)
%   Example: L  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0,'jacobian',JFCN)
%   'disp','none' - sets disply style, 'none' shows no plot (default),
%   'disp','2d' - shows 2d plots time vs values of Lyapunov exponents,
%   'disp','3d' - shows 3d plots with colors corresponding to local values of
%   Lyapunov exponents,
%   'disp','all' - shows all plots
%   Example: L  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0,'disp','3d')
%   'df',N - divides each step by N to calculate the local Lyapunov
%   exponents (default N = 30) 
%   Example: L  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0,'df',10)
%   Output:
%   L - vector of averaged Lyapunov exponents (base e),
%   LSPAN - matrix of local Lyapunov exponents evolution over times TSPAN
%   LEXP - matrix of global Lyapunov exponents evolution over times TSPAN
%-----------------------------------------------------------------------------
% Copyright (C) 2022, Karimov A.I.


fsys = varargin{1,1}; %ODE system derivative function
t = varargin{1,2}; %set time array
y0 = varargin{1,3};

DF = 30; %division factor
fdisp2d = 0; %display 2d
fdisp3d = 0; %display 3d

deltaX = 1e-8; %delta for Jacobian

fsyslin = @(t,D,x)(Jnum(fsys, t, x, deltaX)*D);

for j = 4:2:nargin %in nargin > 3
    if strcmp(varargin{1,j},'jacobian')
        if nargin >= j + 1
            Jfun = varargin{1,j + 1};
            fsyslin = @(t,D,x)(Jfun(t,x)*D);
        else
            error('No agrument matching the property jacobian');
        end
    end
    if strcmp(varargin{1,j},'disp')
        if nargin >= j + 1
            if strcmp(varargin{1,j + 1}, '2d')
                fdisp2d = 1;
            end
            if strcmp(varargin{1,j + 1}, '3d')
                fdisp3d = 1;
            end
            if strcmp(varargin{1,j + 1}, 'all')
                fdisp2d = 1;
                fdisp3d = 1;
            end
        else
            error('No agrument matching the property disp');
        end
    end
    if strcmp(varargin{1,j}, 'df')
        if nargin >= j + 1
            DF = varargin{1,j + 1};
        else
            error('No agrument matching the property df');
        end
    end
end
    
dim =  length(y0);

N = length(t);
Lsum = zeros(dim,1);
lyapexp = zeros(dim,N);

y = zeros(dim,N);
x0 = y0;


del0 = eye(dim); %initial vectors
del1 = del0;

localexp = zeros(dim,N);

xplot = cell(1,N);


hw = waitbar(0,'Calculating Lyapunov Spectrum');

%first, iterate some time before it falls on the attractor
h = t(2) - t(1);


[~,x] = ode45(fsys,0:h:100*h,x0); %fiducial trajectory
x0 = x(end,:);
x0 = transpose(x0);


for i = 2:N
    waitbar(i/N,hw);
    
    h = t(i) - t(i - 1);

    [~,x] = ode45(fsys,t(i - 1):h/DF:t(i),x0); %fiducial trajectory
    x =  transpose(x);
    x1 = x(:,end); %last point
    
    xplot{1,i} = x;
    %for each component in dim
    for k = 1:dim
        [~,dx] = ode45(@(t,d)fsyslin(t,d,x0),t(i - 1):h/DF:t(i),del1(:,k)); %solve linearized system
        del1(:,k) = transpose(dx(end,:));
    end
    
    v = del1;
    
    for k = 1:dim
        vk = v(:,k);
        for l = 1:k-1 %perfrom Gram-Schmidt algorithm step
            u1 = v(:,l);
            vk = vk - (vk'*u1)*u1/(u1'*u1); %subtract projection
        end
        Nk = norm(vk); %compute norm
        localexp(k,i) = 1/h*log(Nk); %local Lyapunov exponent
        v(:,k) = vk/Nk; %normalize k-th component
        
        Lsum(k) = Lsum(k) + localexp(k,i); %sum
        lyapexp(k,i) = Lsum(k)/i; %averaging
    end
    
    del1 = v;
    
    y(:,i) = x1;
    x0  = x1;
end

close(hw);

if fdisp2d

figure(1); hold on
plot(t,y); %draw plot
xlabel('t'),ylabel('Y');

figure(3);
plot(t,lyapexp,t,localexp); %draw plot
xlabel('t'),ylabel('Laypunov exponent');
strleg = cell(1,2*dim);
for k = 1:dim
    strleg{1,k} = ['$\lambda_',num2str(k),' \to $ ',num2str(lyapexp(k,end))];
end
for k = 1:dim
    strleg{1,dim + k} = ['$\lambda_',num2str(k),'$ local'];
end
legend(strleg,'interpreter','latex');
set(gcf,'position',[100    150    1150    350]);
end

% see 3d figures
if fdisp3d
    col1 = [1 0 0]; %red color
    col2 = [0 0 1]; %blue color
    %determine colors
    lmin = zeros(dim,1);
    lmax = zeros(dim,1);
    ldel = zeros(dim,1);
    for k = 1:dim
        lmin(k) = quantile(localexp(k,:),0.15);
        lmax(k) = quantile(localexp(k,:),0.85);
        ldel(k) = lmax(k) - lmin(k);
    end
    
    % initialize figures
    figure(4);
    
    for k = 1:dim
        ax = subplot(1,dim,k); %plot it
        
        Xmat = nan(N,DF+1);
        Ymat = nan(N,DF+1);
        Zmat = nan(N,DF+1);

        colMat = zeros(N,3);
        
        M = 10;

        for i = M:N %do not draw the first M segments of line (inappropriate LLE)
            lk = min([lmax(k),max([localexp(k,i),lmin(k)])]); %local k-th value
            colfactor1 = (lk - lmin(k))/ldel(k);
            colfactor2 = 1 - colfactor1;
            
            x =  xplot{1,i};% extract i-th segment of attractor

            Xmat(i,:) = x(1,:);
            Ymat(i,:) = x(2,:);
            Zmat(i,:) = x(3,:);

            colMat(i,:) = col1*colfactor1 + col2*colfactor2;

        end

        plot3(ax,Xmat',Ymat',Zmat','LineWidth',2); %draw plot
        colororder(ax,colMat);

        %now show colorbar
        nvals = 50;
        cmap = zeros(nvals,3);
        for ctr = 1:nvals
            colfactor1 = (ctr - 1)/(nvals - 1);
            colfactor2 = 1 - colfactor1;
            cmap(ctr,:) = col1*colfactor1 + col2*colfactor2;
        end
        colormap(ax,cmap); %apply colormap
        caxis(ax,[lmin(k) lmax(k)]); %set limits
        colorbar(ax,'Ticks',linspace(lmin(k),lmax(k),5)); %show colorbar with 5 marks
        

        
        view([0 -1 0]);
        xlabel('x'),ylabel('y');zlabel('z');
        title(['$\lambda_',num2str(k),'$'],'interpreter','latex');
        set(gcf,'position',[100    100    1150    350]);

        drawnow;
    end
end

L = lyapexp(:,end); %Estimated spectrum
Lspan = localexp;
Lexp = lyapexp;
end

function J = Jnum(fun, t0, x0, deltaX)
%Jnum finds partial derivatives of order 2

dim = length(x0);

dX = zeros(dim,1);
dX(1) = deltaX;
J = zeros(dim);
for i = 1:dim
    J(:,i) = 0.5*(feval(fun,t0,x0 + dX) - feval(fun,t0,x0 - dX))/deltaX; %ord 2
    dX = circshift(dX,1);
end

end

