function J = Jlorenz2(t,X)
	sig = 10; bet = 8/3; rho = 28;
    x = X(1); y = X(2); z = X(3);
    
    J = [-sig, sig, 0;     % sig*(y - x);
         (rho - z),  -1, -x;     % x*(rho - z) - y;
         y,     x,  -bet]; %x*y - bet*z;

     
end