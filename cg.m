% by Nathaniel Mathews
% A a matrix, b a solution vector and x0 an initial vector.
% Returns a vector x such that Ax = b +/- EPSILON and the number of 
% iterations required to reach that solution.

function [x, iter] = cg(A, b, x0)
    EPSILON = 0.00000001;
    x=x0;
    r=b-A*x;
    p=r;
    rsold=r'*r;
    iter=0;
 
    while r'*r > EPSILON
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        p=r+rsnew/rsold*p;
        rsold=rsnew;
        iter = iter+1;
    end
end