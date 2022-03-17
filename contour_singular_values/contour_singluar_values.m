% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 3/15/22
%
% This script picks two random orthogonal matrices and constructs a new
% matrix with determined singular values by constructing its SVD
% decomposition

% create two random orthogonal matrices
R1 = rand(2);
R2 = rand(2);
[U, ~] = qr(R1);
[V, ~] = qr(R2);

% create A
s1 = 2;
s2 = 1;
A1 = U*[s1 0;0 s2]*V';

% this function has a unique minimum at the origin and will have level
% curves that qualitatively differ as the condition number of A varies
v = [.5; .5];
f1 = @(x, y) [x-v(1); y-v(2)]'*(A1'*A1)*[x-v(1); y-v(2)];

% plot the contours
[X,Y] = meshgrid(-10:.1:10);
z = arrayfun(f1, X, Y);
contour(X,Y,z)
grid on
title('Level Curves of $f(x)=(x-v)^TA^TA(x-v)$, $\frac{\sigma_1}{\sigma_2}=2$', interpreter='Latex')
colorbar

% create another matrix with a larger condition number
figure
s1 = 10;
s2 = .1;
A2 = U*[s1 0;0 s2]*V';
f2 = @(x, y) [x-v(1); y-v(2)]'*(A2'*A2)*[x-v(1); y-v(2)];

% plot the contours
[X,Y] = meshgrid(-10:1:10);
z = arrayfun(f2, X, Y);
contour(X,Y,z)
grid on
title('Level Curves of $f(x)=(x-v)^TA^TA(x-v)$, $\frac{\sigma_1}{\sigma_2}=100$', interpreter='Latex')
colorbar

