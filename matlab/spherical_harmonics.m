% Pass f to data to the function as an MxN array where
% f(i,j,3) = f(theta,phi), theta = f(i,j,1)
% and phi = f(i,j,2)
function coefficients = spherical_harmonics(THETA,PHI,F,n)
% Compute coefficients with spectral method
[M,N] = size(THETA);
coefficients = zeros(n*n,1);
integrand = zeros(M,N);
for l = 0:(n-1)
    for m = -l:l
        for i = 1:M
            for j = 1:N
                theta = THETA(i,j);
                phi = PHI(i,j);
                integrand(i,j) = sin(theta) * F(theta,phi) * conj(harmonicY(l, m, theta, phi,'type','real'));
            end
        end
        % Store coefficients compactly
        idx = int32(l*(l+1) + m + 1);
        coefficients(idx) = trapz(PHI(1,:),trapz(THETA(:,1),integrand,1));
    end
end

end
