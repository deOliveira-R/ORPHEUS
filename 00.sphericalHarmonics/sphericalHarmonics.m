function sphericalHarmonics(l, m)

% Calculates and plots real spherical harmonic functions

num = 50;
theta = 0 : pi/num : pi;
phi = 0 : 2*pi/(num-1) : 2*pi;
for i = 1:num
    for j = 1:num
        P = legendre(l,cos(theta(i)));
        R = sqrt((1*(m==0) + 2*(m~=0))*factorial(l-abs(m))/factorial(l+abs(m))) * ...
            P(abs(m)+1) * ( (m>=0) * cos(m*phi(j)) + (m<0) * sin(abs(m)*phi(j)));

        x(i,j) = abs(R)*sin(theta(i))*cos(phi(j));
        y(i,j) = abs(R)*sin(theta(i))*sin(phi(j));
        z(i,j) = abs(R)*cos(theta(i));
    end
end
f = figure('visible','off');
surf(x,y,z);
axis equal;
axis off;
saveas(f, 'Fig.pdf');
end