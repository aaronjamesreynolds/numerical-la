close all; clear all;
n = 10; a = rand(n,n); q = zeros(n,n); 
r = zeros(n,n); v = zeros(n,n);

for i = 1:n

    v(:,i) = a(:, i);
    
end

for i = 1:n
   
    r(i,i) = norm(v(:,i));
    q(:,i) = v(:,i)/r(i,i);
    
    for j = i + 1:n
        
        r(i, j) = transpose(q(:,i)) * v(:, j);
        v(:,j) = v(:,j) - r(i, j)*q(:,i);
        
    end
 
end
