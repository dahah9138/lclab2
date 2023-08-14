% Code originally written by Benny
% Modified by Darian Hall to be faster
function [chi,ill_defined] = f_chi(nn)

method = 'tensor';

% Eigenvalue error tolerance
tol = 0.01;

[m,n,p,~] = size(nn);

dx = 2/m;
dy = 2/n;
dz = 2/p;

% better than seeing 0, 1, 2,..., m-2 in console
f = waitbar(0, '1', 'Name', 'Generating chi field...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

nx000 = nn(3:end-2,3:end-2,3:end-2,1);
ny000 = nn(3:end-2,3:end-2,3:end-2,2);
nz000 = nn(3:end-2,3:end-2,3:end-2,3);
nx100 = derivativeX(nn(:,:,:,1),dx);
nx010 = derivativeY(nn(:,:,:,1),dy);
nx001 = derivativeZ(nn(:,:,:,1),dz);
ny100 = derivativeX(nn(:,:,:,2),dx);
ny010 = derivativeY(nn(:,:,:,2),dy);
ny001 = derivativeZ(nn(:,:,:,2),dz);
nz100 = derivativeX(nn(:,:,:,3),dx);
nz010 = derivativeY(nn(:,:,:,3),dy);
nz001 = derivativeZ(nn(:,:,:,3),dz);

twist = zeros(m-2,n-2,p-2,3);

x11 = nz100.*ny000 - ny100.*nz000;
x12 = nx100.*nz000 - nz100.*nx000;
x13 = ny100.*nx000 - nx100.*ny000;
x21 = nz010.*ny000 - ny010.*nz000;
x22 = nx010.*nz000 - nz010.*nx000;
x23 = ny010.*nx000 - nx010.*ny000;
x31 = nz001.*ny000 - ny001.*nz000;
x32 = nx001.*nz000 - nz001.*nx000;
x33 = ny001.*nx000 - nx001.*ny000;
ill_defined = zeros(m-4, n-4, p-4);

clear nx000 ny000 nz000 nx100 nx010 ny100 ny010 nz100 nz010 nz001

for i = 1:m-4
    
    simple_waitbar(i, m-4, f)
    
    for j = 1:n-4

     for k = 1:p-4


         Chi_tensor = [x11(i,j,k) x12(i,j,k) x13(i,j,k);...
             x21(i,j,k) x22(i,j,k) x23(i,j,k);...
             x31(i,j,k) x32(i,j,k) x33(i,j,k)]; 

        % Matrix Equation: Chi_tensor * V = V * D

        % V is the matrix of right eigenvectors
        % D is the diagonal eigenvalue matrix
        % W is the matrix of left eigenvectors

        trA2 = trace(Chi_tensor*Chi_tensor);
        tr2A = trace(Chi_tensor)^2;
        discriminant = 2 * trA2 - tr2A;
        
        % Well-defined chi field
        if discriminant > tol
            [V,D] = eig(Chi_tensor);
            [~,eig_label] = max(max(abs(D)));
            twist_slice = V(:,eig_label)';
            
            ill_defined(i,j,k) = 0;
        else
            % Find the circulation direction of the defect line (eigval=0)
            x = null(Chi_tensor);
            ill_defined(i,j,k) = 1;
            twist_slice = x;
        end

        % Normalize
        chi_len = sqrt(dot(twist_slice,twist_slice,1));
        if chi_len == 0
            disp('Faulty tolerance selected')
        end
        twist_slice = twist_slice / chi_len;
        twist(i,j,k,:) = twist_slice(:);

     end
    end
end


f.Name = 'Generated chi field';
delete(f);

chi = real(twist);

% % normalize
modnn = sqrt(chi(:,:,:,1).^2+chi(:,:,:,2).^2+chi(:,:,:,3).^2);
chi_x = chi(:,:,:,1);
chi_y = chi(:,:,:,2);
chi_z = chi(:,:,:,3);
chi_x(modnn>0) = chi_x(modnn>0)./modnn(modnn>0);
chi_y(modnn>0) = chi_y(modnn>0)./modnn(modnn>0);
chi_z(modnn>0) = chi_z(modnn>0)./modnn(modnn>0);

chi(:,:,:,1) = chi_x;
chi(:,:,:,2) = chi_y;
chi(:,:,:,3) = chi_z;

end

function simple_waitbar(ii, sz, f)

waitbar(ii/sz, f, sprintf('%.1f %%', ii/sz*100))

end

function DxF = derivativeX(F, dx)
    % 4th order FD coefficients
    c_0 = 1/12/dx;
    c_1 = -2/3/dx;
    c_2 = 2/3/dx;
    c_3 = -1/12/dx;
    
    DxF = c_3 * F(5:end,3:end-2,3:end-2) + c_2 * F(4:end-1,3:end-2,3:end-2) +...
        c_1 * F(2:end-3,3:end-2,3:end-2) + c_0 * F(1:end-4,3:end-2,3:end-2);
end

function DyF = derivativeY(F, dy)
    % 4th order FD coefficients
    c_0 = 1/12/dy;
    c_1 = -2/3/dy;
    c_2 = 2/3/dy;
    c_3 = -1/12/dy;
    
    DyF = c_3 * F(3:end-2,5:end,3:end-2) + c_2 * F(3:end-2,4:end-1,3:end-2) +...
        c_1 * F(3:end-2,2:end-3,3:end-2) + c_0 * F(3:end-2,1:end-4,3:end-2);
end

function DzF = derivativeZ(F, dz)
    % 4th order FD coefficients
    c_0 = 1/12/dz;
    c_1 = -2/3/dz;
    c_2 = 2/3/dz;
    c_3 = -1/12/dz;
    
    DzF = c_3 * F(3:end-2,3:end-2,5:end) + c_2 * F(3:end-2,3:end-2,4:end-1) +...
        c_1 * F(3:end-2,3:end-2,2:end-3) + c_0 * F(3:end-2,3:end-2,1:end-4);
end
