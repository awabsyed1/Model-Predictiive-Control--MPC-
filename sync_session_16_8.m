%% ACS6116 Synchronous Session 16

%% Problem setup
clear

% System matrices
A = [1 1; 0 2];
B = [0; 0.5];
C = [1 0];

% dimensions
n = size(A,1);
m = size(B,2);

% initial state
x0 = [3; 0];

% horizon N
N = 3;

% cost matrices
Q = C'*C;
R = 1;
K = [-2 -6];
P = dlyap((A+B*K)',Q+K'*R*K);


%% Ex 5.1
umin = -1.5;
umax = +1.5;
xmin = [-10; -5];
xmax = [+10; +5];

% build constraint matrices
Pu = [eye(m); -eye(m)];
qu = [umax; -umin];
Px = [eye(n); -eye(n)];
qx = [xmax; -xmin];
PxN = Px;
qxN = qx;

% build MPC problem matrices
[F, G] = predict_mats(A,B,N);
[H, L, M] = cost_mats(F,G,Q,R,P);
[Pc,qc,Sc] = constraint_mats(F,G,Pu,qu,Px,qx,PxN,qxN);
%[Pc,qc,Sc] = constraint_mats(F,G,Pu,qu,[],[],[],[]);

%% Ex 5.2
nk = 20;
x = x0;
xs = zeros(n,nk);
us = zeros(m,nk);

for k = 1:nk+1
    
    % store x
    xs(:,k) = x;
    
    % solve the MPC problem
    [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);
        
    if flag < 1
        disp(['Optimization infeasible at k = ' num2str(k)])
    break 
    end
    
    % extract the first control
    u = Uopt(1:m);
    
    % store u
    us(:,k) = u;
    
    % apply to the plant
    x = A*x + B*u;
    
end

figure(1)
plot((0:k-1),xs(:,1:k)')

figure(2)
stairs((0:k-1),us(:,1:k)')

%% Ex 5.2(h)

PxN = Px;
qxN = 0*qx;

% build MPC problem matrices
[F, G] = predict_mats(A,B,N);
[H, L, M] = cost_mats(F,G,Q,R,P);
[Pc,qc,Sc] = constraint_mats(F,G,Pu,qu,Px,qx,PxN,qxN);

figure(3)
hold on

for x1 = -10:0.2:10
    for x2 = -5:0.1:5
        
        x = [x1; x2];
        
        % solve the MPC problem
        [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);
        
        if flag >= 1
            plot(x1,x2,'b.');
        end
        
    end
end
