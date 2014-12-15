% by Nathaniel Mathews
% A a matrix.
% Returns Q an orthogonal matrix and R an
% upper triangular matrix such that A = QR.
% Uses the Householder Reflectors method.

function [Q, R] = qr(A)
    if max(size(A)) == 1
        R = A
        Q = 1
        return
    end
    I = zeros(size(A,1), size(A,1));
    for n=[1:min(size(I))]
        I(n,n) = 1;
    end
    x = A(:,1);
    w = zeros(size(x));
    w(1) = norm(x);
    u = x - w;
    P = (u*transpose(u))/(transpose(u)*u)
    H = I - 2*P
    A = H*A
    Htilde = 1;
    if min(size(A)) >= 2
        [Htilde, A(2:size(A,1),2:size(A,2))] =...
            qr(A(2:size(A,1),2:size(A,2)))
    end
    Q = zeros(size(H));
    Q(1,1) = 1;
    Q(2:size(Q,1), 2:size(Q,2)) = Htilde;
    Q = H*Q
    R = A
    return
end