function[Q, R] = mgs(A)
    
    n = size(A,2);

    Q = zeros(size(A)); 
    R = zeros(n); v = zeros(size(A));
        
    for i = 1:n
        v(:,i) = A(:, i);
    end

    for i = 1:n
        R(i,i) = norm(v(:,i));
        Q(:,i) = v(:,i)/R(i,i);

        for j = i + 1:n
            R(i, j) = transpose(Q(:,i)) * v(:, j);
            v(:,j) = v(:,j) - R(i, j)*Q(:,i);
        end
    end
end
