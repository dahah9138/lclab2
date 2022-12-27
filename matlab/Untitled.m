x = -3:.5:3; 
y = -3:.5:3; 
[X,Y] = meshgrid(x,y);
F = X.^2 + Y.^2;
I = trapz(y,trapz(x,F,2))