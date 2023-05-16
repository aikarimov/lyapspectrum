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
%   'trans',TTRANS skips TTRANS time before calculating the spectrum
%   'df',N - divides each step by N to calculate the local Lyapunov
%   exponents (default N = 30) 
%   Example: L  = LYAPSPECTRUM(ODEFUN,TSPAN,Y0,'df',10)
%   Output:
%   L - vector of averaged Lyapunov exponents (base e),
%   LSPAN - matrix of local Lyapunov exponents evolution over times TSPAN
%   LEXP - matrix of global Lyapunov exponents evolution over times TSPAN
%-----------------------------------------------------------------------------
% Copyright (C) 2023, Karimov A.I.


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

y(:,1) = x0;

%main cycle
h = t(2) - t(1); %suppose, stepsize is uniform
[~,xfid] = ode78(fsys,0:h/DF:t(end) + Ttrans,x0); %fiducial trajectory

for i = 2:N
    waitbar(i/N,hw);
    
    %h = t(i) - t(i - 1);
    ilbnd = (i-2)*DF + 1;
    iubnd = (i-1)*DF + 1;
    x = xfid(ilbnd:iubnd,:);
    x =  transpose(x);
    x1 = x(:,end); %last point
    
    xplot{1,i} = x;
    %for each component in dim
    for k = 1:dim
        [~,dx] = ode78(@(t,d)fsyslin(t,d,x0),t(i - 1):h/DF:t(i),del1(:,k)); %solve linearized system
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


%cut off transient
Nt = max([ceil(Ttrans/h),1]);
y = y(:,Nt:end);
lyapexp = lyapexp(:,Nt:end);
localexp = localexp(:,Nt:end);
t = t(Nt:end) - t(Nt);

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

%3d figures - new version

if fdisp3d
    if dim < 3
        disp('Dimension is less then 3');
    else
        FontSIZE = 10;
        figure(figctr + 2);
        for k = 1:dim
            subplottight(1,dim,k);
            colormap(turbo)
            surf([y(1,:); y(1,:)], [y(2,:); y(2,:)], [y(3,:); y(3,:)],[localexp(k,:);localexp(k,:)],'EdgeColor','flat', 'FaceColor','none')
            set(gca,'FontSize',FontSIZE-1);
            title(['$\lambda_',num2str(k),'$'],'interpreter','latex');
            set(gca, 'XDir', 'reverse'); hold on;
            zlabel('$z$','interpreter','latex','FontSize',FontSIZE);
            ylabel('$y$','interpreter','latex','FontSize',FontSIZE);
            xlabel('$x$','interpreter','latex','FontSize',FontSIZE);
            set(gca,'TickLabelInterpreter','latex');
            c = colorbar;
            set(c,'TickLabelInterpreter','latex');
            view(viewvect);
        end
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

function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    ax = subplot('Position', [(c-0.82)/m, 1-(r-0.2)/n, (0.75)/m, (0.7)/n]);
    if(nargout > 0)
      h = ax;
    end
end