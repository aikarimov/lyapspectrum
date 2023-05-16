function J = Jjosep(t,x)
%JJOSEP - Jacobian for ODE of the shunted Josephson junction
bL = 29.215;
bC = 0.707;
VgIcRs = 6.9;
RsRN = 0.367;
RsRsg = 0.0478;

if (abs(x(2)) > VgIcRs)
    g = RsRN;
else
    g = RsRsg;
end

J = zeros(3);
J(1,:) = [0, 1, 0];
J(2,:) = 1/bC*[-cos(x(1)),  -g, -1];   %1/bC*(i - g*x(2) - sin(x(1)) - x(3));
J(3,:) = [0, 1/bL, -1/bL];

end