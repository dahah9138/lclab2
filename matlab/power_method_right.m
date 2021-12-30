% Darian Hall
% Computes the dominant right eigenvector and eigenvalue of M with specified tolerance
function [eigenvector, eigenvalue] = power_method_right(M, tolerance)

    % Initial Guess
    eigenvector = 0.5*ones(size(M(:,1)));
    eigenvector(1) = 1;
    
    eigenvalue = 1;
    
    loop = 1;
    
    while(loop)
        
        % Remember previous eigenvalue
        prev_eigenvalue = eigenvalue;
        
        % Compute new eigenvector and eigenvalue
        eigenvector = M * eigenvector;
        eigenvalue = max(eigenvector);
        eigenvector = eigenvector / eigenvalue;
        
        % Check if we should continue
        loop = (abs(prev_eigenvalue - eigenvalue) > tolerance);
        
    end
    
end