% system dimensions
n = 3; % no. of states
m = 1; % no. of inputs

% input limits
umax = 1;
umin = -2;

% state limits
xmax = [2; 1; 4];
xmin = [-3; -2; -5];

% matrices Pu*u <= qu
% 1*u <= 1 (row 1)
% -1*u <= 2 (row 2) [which implies 1*u >= -2]
Pu = [eye(m); -eye(m)];
qu = [umax; -umin];
Px = [eye(n); -eye(n)];
qx = [xmax; -xmin];
Pxf = Px;
qxf = qx;

plotregion(-Px,-qx)