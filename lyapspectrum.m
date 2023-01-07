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
%   'view',VECT sets 3d display with the view defined by VECT ([0 0 1] by
%   default)
%   'transient',TTRANS skips TTRANS time before calculating the spectrum
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

viewvect = [0 0 1];
Ttrans = 0;

deltaX = 1e-8; %delta for Jacobian

fsyslin = @(t,D,x)(Jnum(fsys, t, x, deltaX)*D); %if we have no analytical Jacobian

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
    if strcmp(varargin{1,j}, 'view')
        if nargin >= j + 1
            viewvect = varargin{1,j + 1};
        else
            error('No agrument matching the property view');
        end
    end
     if strcmp(varargin{1,j}, 'trans')
        if nargin >= j + 1
            Ttrans = varargin{1,j + 1};
        else
            error('No agrument matching the property trans');
        end
    end
end

figctr = 1; %figure counter

if(fdisp3d || fdisp2d) %if any figure is shown
     hw =  findobj('type','figure'); %find all figures
     hwctr = length(hw);
     Nmax = 0;
     for j = 1:hwctr %find maximal figure number
        numb = hw(j,1).Number;
        if numb > Nmax
            Nmax = numb;
        end
     end
     figctr = Nmax + 1;
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
if Ttrans > 0
    h = t(2) - t(1);
    [~,x] = ode45(fsys,0:h:Ttrans,x0); x = x'; %fiducial trajectory
    x0 = x(end,:);
    x0 = transpose(x0);
end

%main cycle
h = t(2) - t(1); %suppose, stepsize is uniform
[~,xfid] = ode78(fsys,0:h/DF:t(end),x0); %fiducial trajectory
%[~,xfid] = DOPRI78s(fsys,0:h/DF:t(end),x0); xfid = xfid';%fiducial trajectory

for i = 2:N
    waitbar(i/N,hw);
    
    %h = t(i) - t(i - 1);
    %[~,x] = ode45(fsys,t(i - 1):h/DF:t(i),x0); %fiducial trajectory
    ilbnd = (i-2)*DF + 1;
    iubnd = (i-1)*DF + 1;
    x = xfid(ilbnd:iubnd,:);
    x =  transpose(x);
    x1 = x(:,end); %last point
    
    xplot{1,i} = x;
    %for each component in dim
    for k = 1:dim
        [~,dx] = ode78(@(t,d)fsyslin(t,d,x0),t(i - 1):h/DF:t(i),del1(:,k)); %solve linearized system
        %[~,dx] = DOPRI78s(@(t,d)fsyslin(t,d,x0),t(i - 1):h/DF:t(i),del1(:,k)); dx = dx';%solve linearized system
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

figure(figctr); hold on
plot(t,y); %draw plot
xlabel('t'),ylabel('Y');

figure(figctr + 1);
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
    col3 = [0 1 0]; %green color
    %determine colors
    lmin = zeros(dim,1);
    lmax = zeros(dim,1);
    ldel = zeros(dim,1);
    for k = 1:dim %for each dimension, set own limits
        lmin(k) = quantile(localexp(k,:),0.15); %below lmin, all is blue
        lmax(k) = quantile(localexp(k,:),0.85); %under lmax, all is red
        ldel(k) = lmax(k) - lmin(k); %interval
    end
    
    % initialize figures
    figure(figctr + 2);
    
    for k = 1:dim
        ax = subplot(1,dim,k); %plot it
        
        Xmat = nan(N,DF+1);
        Ymat = nan(N,DF+1);
        Zmat = nan(N,DF+1);

        colMat = zeros(N,3);
        
        M = 10;

        for i = M:N %do not draw the first M segments of line (inappropriate LLE)
            lk = min([lmax(k),max([localexp(k,i),lmin(k)])]); %local k-th value     
            x =  xplot{1,i};% extract i-th segment of attractor

            Xmat(i,:) = x(1,:);
            Ymat(i,:) = x(2,:);

            if dim > 2 %if dimensions are more than 2, we set variable Z
                Zmat(i,:) = x(3,:); 
            end

            colMat(i,:) = assigncolor(col1,col2,col3,lmin(k),lmax(k),ldel(k),lk);

        end
        if dim > 2
            plot3(ax,Xmat',Ymat',Zmat','LineWidth',2); %draw plot
        else
            plot(ax,Xmat',Ymat','LineWidth',2); %draw plot
        end
        colororder(ax,colMat);

        %now show colorbar
        nvals = 50; %values of color
        cmap = zeros(nvals,3);
        lvals = linspace(lmin(k),lmax(k),nvals);
        for ctr = 1:nvals
            cmap(ctr,:) = assigncolor(col1,col2,col3,lmin(k),lmax(k),ldel(k),lvals(ctr));
        end
        colormap(ax,cmap); %apply colormap
        caxis(ax,[lmin(k) lmax(k)]); %set limits
        colorbar(ax,'Ticks',linspace(lmin(k),lmax(k),5),'TickLabelInterpreter','latex'); %show colorbar with 5 marks
        view(viewvect);
        xlabel('$x_{1}$','interpreter','latex'),ylabel('$x_{2}$','interpreter','latex');
        if dim >= 3
            zlabel('$x_{3}$','interpreter','latex');
        end
        set(ax,'TickLabelInterpreter','latex');
        title(['$\lambda_',num2str(k),'$'],'interpreter','latex');
        set(gcf,'position',[100    100    1350    350]);

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

function col = assigncolor(col1,col2,colmid,lmin,lmax,ldel,lk)
%col1 for max
%col2 for min
%colmid for middle value
%lmin min value
%lmax max value
%ldel interval

if lmin*lmax < 0 %if signs are different
    ldel2 = -lmin; %set middle color to zero values
    ldiv = 0;
else
    ldel2 = ldel/2; %set middle color to exact middle
    ldiv = lmin + ldel2;
end
if lk > ldiv
    colfactor1 = (lk - ldiv)/(lmax - ldiv);
    colfactor2 = 1 - colfactor1;
    %assign color
    col = col1*colfactor1 + colmid*colfactor2;
else
    colfactor1 = (lk - lmin)/ldel2;
    colfactor2 = 1 - colfactor1;
    %assign color
    col = colmid*colfactor1 + col2*colfactor2;
end
%set limits
for i = 1:3
    if col(i) < 0 
        col(i) = 0;
    end
    if col(i) > 1
        col(i) = 1;
    end
end
end