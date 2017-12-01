function x = mysolve(A, b, m, flag)
% Fast solution of A*x = b where A is banded + m dense rows via Schur
% complement factorisation.

n = size(A);

if ( n <= m )
    x = A\b;
    return
end

i1 = 1:m;
i2 = m+1:n;
i3 = 2:m+1;

if ( nargin < 4 || flag == 0 || isempty(i2) )
    c = A(i2,i2)\[b(i2), A(i2,i1)];
else
    % Row scaling to improve accuracy
    AA = A(i2,i2);
    s = 1./ max(1, max(abs(AA), [], 2) );  
    AA = bsxfun(@times, s, AA);
    bb = s.*[b(i2), A(i2,i1)];
    c = AA\bb;
end

x = (A(i1,i1) - A(i1,i2)*c(:,i3)) \ (b(i1) - A(i1,i2)*c(:,1));
x = [x ; c*[1; -x]];
end