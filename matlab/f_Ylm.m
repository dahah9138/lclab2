% Function that takes f_lm coefficients of a function f(theta,phi)
function f = f_Ylm(theta, phi, coeffs, n)
    % Evaluate the function
    f = 0;
    for l=0:(n-1)
        for m=-l:l
            id = l*(l+1)+m+1;
            f=f+coeffs(id)*harmonicY(l,m,theta,phi,'type','real');
        end
    end
end