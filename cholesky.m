% by Nathaniel Mathews
% A a Hermitian, Positive-Definite matrix.
% Returns a lower triangular matrix L such that
% LL* = A.

function L = cholesky(A)
    L = zeros(size(A));
    for i=[1:size(A,1)]
        L(i,i) = 1;
    end
    
    for i=[1:size(A,1)-1]
        Li = zeros(size(A));
        Li((i+1):size(A,1), i) = A((i+1):size(A,1), i) / sqrt(A(i,i));
        
        Anext = zeros(size(A));
        
        for j=[1:size(A,1)]
            if j == i
                Li(j,j) = sqrt(A(i,i));
            else
                Li(j,j) = 1;
            end
            Anext(j,j) = 1;
        end
        
        Anext((i+1):size(A,1), (i+1):size(A)(2)) =...
            A((i+1):size(A,1), (i+1):size(A,2)) - 1 / A(i,i)...
            * A((i+1):size(A,1),i) * A(i, (i+1):size(A,2));
        Li = Li
        L = L*Li
        A = Anext
    end
    L(size(A,1),size(A,2)) = sqrt(A(size(A,1),size(A)(2)));
    return
end