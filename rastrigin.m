%% rastrigin's function
x=-6:0.01:6;
y=-6:0.01:6;
[X,Y]=meshgrid(x,y);
F=20+X.^2+Y.^2-10*(2*cos(2*pi*X)+2*cos(2*pi*Y));
s=surf(X,Y,F)
s.EdgeColor='none'