% by Nathaniel Mathews
% A a matrix.
% Returns Q an orthogonal matrix and R an
% upper triangular matrix such that A = QR.

function [Q,R] = classical_gramschmidt_qr(A)
    Q = zeros(size(A));
    R = zeros(size(A,2));
    R(1,1) = (transpose(A(:,1)) * A(:,1))^(1/2);
    Q(:,1) = A(:,1) / R(1,1);
    
    for j=[2:size(A,2)]
        y = A(:,j);
        for i=[1:(j-1)]
            R(i,j) = transpose(A(:,j)) * Q(:,i);
            y = y - R(i,j) * Q(:,i);
        end
        R(j,j) = (transpose(y) * y) ^ (1/2);
        Q(:,j) = y / R(j,j);
    end
    
    return
end