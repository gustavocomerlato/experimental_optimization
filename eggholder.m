%% Eggholder function
x=-512:1:512;
y=-512:1:512;
[X,Y]=meshgrid(x,y);
F=-(Y+47).*sin(sqrt(abs(X./2+Y+47)))-X.*sin(sqrt(abs(X-(Y+47))));
s=surf(X,Y,F)
s.EdgeColor='none'
figure(2)
contourf(X,Y,F)