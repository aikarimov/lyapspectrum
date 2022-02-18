function dX = lorenz2(t,X)
    sig = 10; bet = 8/3; rho = 28;
    dX = X;
    x = X(1); y = X(2); z = X(3);
    dX(1) = sig*(y - x);
    dX(2) = x*(rho - z) - y;
    dX(3) = x*y - bet*z;
end