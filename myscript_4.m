%% Problem setup

% system matrices
A = [3 1; 2.01 1.99];
B = [0.1; 2.1];
C = [-0.35 1];

% controller parameters
Q = C'*C;

% change these
R = 100;
P = 0*Q;
N = 3;

%% Synchronous session 6
% prediction matrices
[F,G] = predict_mats(A,B,N);

% cost matrices
 [H,L,M] = cost_mats(F,G,Q,R,P);

% optimal policy
S = -H\L;

% optimal feedback law
KN = S(1,:);

% form closed-loop system
Phi = A+B*KN;

% stability check
rho = max(abs(eig(Phi)));

if rho >= 1
    display('System is not stable')
else
    display('System is stable')
end

%% Synchronous session 7

% dimensions
n = size(A,2);
m = size(B,2);

% constraints
Px = zeros(n,n*N);
Px(:,(end-n+1):end) = eye(n);

% E and D, so that E*U = -D*x
E = Px*G;
D = Px*F;

% KKT system (LHS)
curlyK = [H E'; E zeros(n)];

% solve KKT system curlyK*[sol] = -[L; D]*x
T = -curlyK\[L; D];

% so that solution [uopt; lambdaopt] = T*x
KNbar = T(1:m,:);

% form closed-loop system
Phibar = A+B*KNbar;

% stability check
rho = max(abs(eig(Phibar)));

if rho >= 1
    display('System with terminal constraint is not stable')
else
    display('System with terminal constraint is stable')
end

%%  Synchronous session 8

% find a stabilizing K for (A,B)
%K = -dlqr(A,B,Q,R);
desired_poles = [0 0];
K = -acker(A,B,desired_poles);

% solve Lyapunov equation
P = dlyap((A+B*K)',Q + K'*R*K);

% cost matrices
[H,L,M] = cost_mats(F,G,Q,R,P);

% optimal policy
S = -H\L;

% optimal feedback law
KN = S(1,:);

% form closed-loop system
Phi = A+B*KN;

% stability check
rho = max(abs(eig(Phi)));

if rho >= 1
    display('System with terminal P is not stable')
else
    display('System with terminal P is stable')
end

%% Synchronous session 10

% initial state
x = [-3; 0];

% Obtain a plot of x2(k+i|k) versus x1(k+i|k)
%Uopt = S*x;
%Xopt = F*x + G*Uopt;
Xopt = (F + G*S)*x;

% Xopt contains x(k+1|k), x(k+2|k), ..., x(k+N|k)
Xopt = [x; Xopt];

% extract x1 and x2 elements
Xoptodd = Xopt(1:2:end,:);
Xopteven = Xopt(2:2:end,:);



% x(k+N|k) is final 2 elements of Xopt
xN = Xopt(end-n+1:end,:);

% x(k+N+1|k) = A*x(k+N|k) + B*K*x(k+N|k) = (A+B*K)*x(k+N|k)
% x(k+N+2|k) = (A+B*K)*x(k+N+1|k) = (A+B*K)^2 * x(k+N|k)
%
% x(k+N+i|k) = (A+B*K)^i * x(k+N|k)

% mode-2 horizon (how many steps I want to see)
N2 = 2;
[F2, G2] = predict_mats(A+B*K,B,N2);

% F2 is [(A+BK); (A+BK)^2; ...; (A+BK)^N2]

% mode 2 state predictions
Xopt2 = [xN; F2*xN;];
% extract x1 and x2 elements
Xopt2odd = Xopt2(1:2:end,:);
Xopt2even = Xopt2(2:2:end,:);

% plot mode 1
figure
plot(Xoptodd,Xopteven,'bo-'), hold on
plot(Xopt2odd,Xopt2even,'r*-')
