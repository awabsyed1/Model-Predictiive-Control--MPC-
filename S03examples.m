%% S03EXAMPLES.M
%
% Script to support ACS6116 synchronous session 3.
%
% This is rough and uncommented code -- use at your own risk! I take no
% responsibility for errors.
%
% (c) P. Trodden, 2021.

%% 1: Q, R ≻ 0 and (A, B) stabilizable.
display('Example 1')
A = [1 1; 0 0.9];
B = [1; 0];
Q = eye(2);
R = 1;

check_ABQR(A,B,Q,R);

%% 2: Q ≽ 0, R ≻ 0, (A, B) reachable, but (Q1/2, A) not detectable
display('Example 2')
A = [1 1; 0 1];
B = [1; 1];
Q = 0*eye(2);
R = 1;

check_ABQR(A,B,Q,R);

%% 3: Q ≽ 0, R ≻ 0, (A, B) stabilizable, and (Q1/2, A) detectable
display('Example 3')
A = [0.9 1; 0 0.9];
B = [1; 0];
Q = 0*eye(2);
R = 1;

check_ABQR(A,B,Q,R);

%% 4: Q, R ≻ 0, (A, B) not stabilizable
display('Example 4')
A = [1 1; 0 1.1];
B = [1; 0];
Q = eye(2);
R = 1;

check_ABQR(A,B,Q,R);

%% 5: Q, R ≽ 0, (A, B) reachable
display('Example 5')
A = [0.9 1; 0 0.9];
B = [0; 1];
Q = 0*eye(2);
R = 0*1;

check_ABQR(A,B,Q,R);

%% 6: Q ≽ 0, R ≻ 0, (A, B) controllable, and (Q1/2, A) detectable
display('Example 6')
A = [1 1; 0 0];
B = [1; 0];
C = [0 1];
Q = C'*C;
R = 1;

check_ABQR(A,B,Q,R);

%% function for checking properties
check_ABQR(A,B,Q,R);

function check_ABQR(A,B,Q,R)

if min(eig(R)) > 0
    display('R is positive definite')
elseif min(eig(R)) == 0
    display('R is positive semidefinite -- DLQR.M will fail!')
else
    display('R is not positive semidefinite  -- DLQR.M will fail!')
end

if min(eig(Q)) > 0
    display('Q is positive definite')
elseif min(eig(Q)) == 0
    display('Q is positive semidefinite')
    % check detectability
    [ABAR,BBAR,CBAR,T,K] = obsvf(A,B,sqrtm(Q));
    nobs = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nobs,1:nobs)))) == 0
        display('(Q^0.5,A) is constructible')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(Q^0.5,A) is detectable')
    else
        display('(Q^0.5,A) is not detectable -- DLQR.M will fail!')
    end
    
    
else
    display('Q is not positive semidefinite -- DQLR.M will fail!')
end

r = rank(ctrb(A,B));
if r == 2
    display('(A,B) is reachable')
else
    display('(A,B) is not reachable')
    % check stabilizability
    [ABAR,BBAR,CBAR,T,K] = ctrbf(A,B,sqrtm(Q));
    nrea = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nrea,1:nrea)))) == 0
        display('(A,B) is controllable')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(A,B) is stabilizable')
    else
        display('(A,B) is not stabilizable -- DLQR.M will fail!')
    end
end

end