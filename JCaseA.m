function J = JCaseA(t,X)
    x = X(1); y = X(2); z = X(3);

%     dX(1) = y;
%     dX(2) = -x + y*z;
%     dX(3) = 1 - y^2;
    
    J = [0, 1, 0; 
         -1, z, y; 
         0,     -2*y,  0]; 
     
end