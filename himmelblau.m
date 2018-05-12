%% Himmelblau's function
x=-5:0.01:5;
y=-5:0.01:5;
[X,Y]=meshgrid(x,y);
F=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
s=surf(X,Y,F)
s.EdgeColor='none'