clear
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

% find a stabilizing K for (A,B)
%K = -dlqr(A,B,Q,R);
desired_poles = [0 0];
K = -acker(A,B,desired_poles);

% solve Lyapunov equation
P = dlyap((A+B*K)',Q + K'*R*K);

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
    display('System with terminal P is not stable')
else
    display('System with terminal P is stable')
end





