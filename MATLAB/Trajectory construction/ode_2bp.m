function dy = ode_2bp(~, y, mu)
dy = zeros(size(y));
r = norm([y(1) y(2) y(3)]);
dy(1:3) = y(4:6); 
dy(4) = (-mu/r^3)*y(1);
dy(5) = (-mu/r^3)*y(2);
dy(6) = (-mu/r^3)*y(3);