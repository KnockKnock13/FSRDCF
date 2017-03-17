function [ f ] = gauss_seidel( x,y,intf,num_GS_iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global WW Y 
%     num_GS_iter = 4;

    
    y = sparse(double(bsxfun(@times, Y, conj(y))));%)
    X = x;
    
    A =  diag(X(:)) + WW;
    AL = tril(A);
    AU = A - AL;
    
%     f1 = intf_r;
%     for iter = 1:num_GS_iter
%         f1= AL \ (b1(:) - AU  * f1);
%     end
%     
%     f2= intf_i;
%     for iter = 1:num_GS_iter
%         f2= AL \ (b2(:) - AU  * f2);
%     end
%     f = real(f1) + 1i*real(f2);
     f = intf;
    for iter = 1:num_GS_iter
        f= AL \ (y(:) - AU  * f);
    end
end

