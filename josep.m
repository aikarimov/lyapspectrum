function dx = josep(t,x)
%JOSEP - ODE of the shunted Josephson junction
bL = 29.215; % 2.6; 
bC = 0.707;
i = 1.25;
VgIcRs = 6.9;
RsRN = 0.367;
RsRsg = 0.0478;

if (abs(x(2)) > VgIcRs)
    g = RsRN;
else
    g = RsRsg;
end

dx = x;

dx(1) = x(2);
dx(2) = 1/bC*(i - g*x(2) - sin(x(1)) - x(3));
dx(3) = 1/bL*(x(2) - x(3));

end