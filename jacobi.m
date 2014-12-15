% by Nathaniel Mathews
% A a matrix, b a solution vector, x0 an initial vector
% returns the vector x such that Ax = b +/- EPSILON, and the number
% of iterations necessary to obtain such a vector.

function [x,iters] = jacobi(A, b, x0)
    EPSILON = 0.00000001;
    x = A*x0;
    iters = 0;
    while transpose((A*x-b)) * (A*x-b) > EPSILON
        U = zeros(size(A));
        L = zeros(size(A));
        Dinv = zeros(size(A));
        for i=[1:min(size(A))]
            Dinv(i,i) = 1/A(i,i);
            U(1:i,i) = A(1:i,i);
            L(i:size(A,1), i) = A(i:size(A,1), i);
            U(i,i) = 0;
            L(i,i) = 0;
        end
        x = Dinv * b - Dinv * (L+U) * x;
        iters = iters+1;
    end
    return
end