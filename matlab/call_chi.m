% Code originally written by Benny
% Modified by Darian Hall to be faster
clear
clc

file_loc = "";
file = "dimer";
file_type = ".mat";

method = 'tensor';

% Eigenvalue error tolerance
tol = 0.0001;

try
    load(strcat(file_loc, file, file_type), 'nn');
catch
    disp("Failed to find nn")
end

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


switch(method)
    case 'tensor'
    x11 = nz100.*ny000 - ny100.*nz000;
    x12 = nx100.*nz000 - nz100.*nx000;
    x13 = ny100.*nx000 - nx100.*ny000;
    x21 = nz010.*ny000 - ny010.*nz000;
    x22 = nx010.*nz000 - nz010.*nx000;
    x23 = ny010.*nx000 - nx010.*ny000;
    x31 = nz001.*ny000 - ny001.*nz000;
    x32 = nx001.*nz000 - nz001.*nx000;
    x33 = ny001.*nx000 - nx001.*ny000;
    eigenvalue = zeros(m-4, n-4, p-4);
    
    clear nx000 ny000 nz000 nx100 nx010 ny100 ny010 nz100 nz010 nz001

    for i = 1:m-4
        
        simple_waitbar(i, m-4, f)
        
        for j = 1:n-4

         for k = 1:p-4


             Chi_tensor = [x11(i,j,k) x12(i,j,k) x13(i,j,k);...
                 x21(i,j,k) x22(i,j,k) x23(i,j,k);...
                 x31(i,j,k) x32(i,j,k) x33(i,j,k)]; 



            % V is the matrix of right eigenvectors
            % D is the diagonal eigenvalue matrix
            % W is the matrix of left eigenvectors

            % Matrix Equation: W' * Chi_tensor = D * W'

            %[~,D,W] = eig(Chi_tensor);
            %[eig_value,eig_label] = max(max(real(D)));

            [eig_vec, eig_value] = power_method_left(Chi_tensor, tol);
            
            eigenvalue(i,j,k) = eig_value;

            %twist_slice = eig_value * W(:,eig_label)';
            twist_slice = eig_value * eig_vec;

            twist(i,j,k,:) = twist_slice(:);

         end
        end
    end
    
    case 'vector'
        % n cross (u dot \nabla)n ~ chi 
        % Rev. Mod. Phys. 46, 617 (1974)
        % not general 
        
        % u = x unit vector
        twistux_x = ny000.*nz100 - nz000.*ny100;
        twistux_y = nz000.*nx100 - nx000.*nz100;
        twistux_z = nx000.*ny100 - ny000.*nx100;
        % u = y unit vector
        twistuy_x = ny000.*nz010 - nz000.*ny010;
        twistuy_y = nz000.*nx010 - nx000.*nz010;
        twistuy_z = nx000.*ny010 - ny000.*nx010;
        % u = z unit vector
        twistuz_x = ny000.*nz001 - nz000.*ny001;
        twistuz_y = nz000.*nx001 - nx000.*nz001;
        twistuz_z = nx000.*ny001 - ny000.*nx001;
        
        modtwstx = sqrt(twistux_x.^2+twistux_y.^2+twistux_z.^2);
        modtwsty = sqrt(twistuy_x.^2+twistuy_y.^2+twistuy_z.^2);
        modtwstz = sqrt(twistuz_x.^2+twistuz_y.^2+twistuz_z.^2);
        
        % compare modtwst and substitute components with the largest one
        % compare u = x, y
        twistux_x(modtwstx<modtwsty) = twistuy_x(modtwstx<modtwsty);
        twistux_y(modtwstx<modtwsty) = twistuy_y(modtwstx<modtwsty);
        twistux_z(modtwstx<modtwsty) = twistuy_z(modtwstx<modtwsty);
        
        % compare u = x, z
        twistux_x(modtwstx<modtwstz) = twistuz_x(modtwstx<modtwstz);
        twistux_y(modtwstx<modtwstz) = twistuz_y(modtwstx<modtwstz);
        twistux_z(modtwstx<modtwstz) = twistuz_z(modtwstx<modtwstz);
        
        twist = cat(4, twistux_x, twistux_y, twistux_z);
%         nn = cat(4, twistuz_x, twistuz_y, twistuz_z);    
%         nn = cat(4, twistuy_x, twistuy_y, twistuy_z);  

        % % normalize
        modnn = sqrt(twist(:,:,:,1).^2+twist(:,:,:,2).^2+twist(:,:,:,3).^2);
        twist(:,:,:,1) = twist(:,:,:,1)./modnn;
        twist(:,:,:,2) = twist(:,:,:,2)./modnn;
        twist(:,:,:,3) = twist(:,:,:,3)./modnn;
end
%%

f.Name = 'Generated chi field';
delete(f);

chi = real(twist);
%chi = abs(twist);



% % normalize
modnn = sqrt(chi(:,:,:,1).^2+chi(:,:,:,2).^2+chi(:,:,:,3).^2);
chi(:,:,:,1) = chi(:,:,:,1)./modnn;
chi(:,:,:,2) = chi(:,:,:,2)./modnn;
chi(:,:,:,3) = chi(:,:,:,3)./modnn;


save(strcat(file_loc, file, '_chi', file_type), 'chi');

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
