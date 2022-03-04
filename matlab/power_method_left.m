% Darian Hall
% Computes the dominant left eigenvector and eigenvalue of M with specified tolerance
function [eigenvector, eigenvalue] = power_method_left(M, tolerance, stop_after)

    % Initial Guess
    eigenvector = 0.5*ones(size(M(1,:)));
    eigenvector(1) = 1;
    
    eigenvalue = 1;
    
    loop = 1;
    
    count = 0;
    
    while(loop)
        
        % Remember previous eigenvalue
        prev_eigenvalue = eigenvalue;
        
        % Compute new eigenvector and eigenvalue
        eigenvector = eigenvector * M;
        eigenvalue = max(eigenvector);
        eigenvector = eigenvector / eigenvalue;
        
        % Check if we should continue
        loop = (abs(prev_eigenvalue - eigenvalue) > tolerance) && (abs(prev_eigenvalue + eigenvalue) > tolerance);
        count = count + 1;
        
        if count > stop_after
            loop = 0;
        end
        
    end
    
end