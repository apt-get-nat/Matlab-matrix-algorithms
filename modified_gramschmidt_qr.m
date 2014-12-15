% by Nathaniel Mathews
% A a matrix.
% Returns Q an orthogonal matrix and R an
% upper triangular matrix such that A = QR.

function [Q,R] = modified_gramschmidt_qr(A)
    Q = A;
    R = zeros(size(A,2));
    for i=[1:size(A,2)]
        R(i,i) = (transpose(Q(:,i)) * Q(:,i))^(1/2);
        Q(:,i) = Q(:,i) / R(i,i);
        for j = [i+1:size(A,2)]
            R(i,j) = transpose(A(:,j)) * Q(:,i);
            Q(:,j) = Q(:,j) - R(i,j) * Q(:,i);
        end
    end
    return
end