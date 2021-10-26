clc;clear;

syms q1 q2 q3 q4 dq1 dq2 dq3 dq4 ddq1 ddq2 ddq3 ddq4 u1 u2 u3
syms m1 m2 m3 m4 r1 r2 L1 L2 J1 J2 g

q = [q1 q2 q4].';
dq = [dq1 dq2 dq4].';
ddq = [ddq1 ddq2 ddq4].';

p1 = [r1*q1, r1];
p2 = p1 + [L1*sin(q2), L1*cos(q2)];
p4 = p2 + [L2*sin(q2), L2*cos(q2)];
v1 = jacobian(p1,q)*dq;
v2 = jacobian(p2,q)*dq;
v4 = jacobian(p4,q)*dq;
KE1 = 1/2*m1*(v1.'*v1) + 1/2*J1*(dq1.'*dq1);
KE2 = 1/2*(m2+m3)*(v2.'*v2);
KE4 = 1/2*m4*(v4.'*v4) + 1/2*J2*((dq4+dq2).'*(dq4+dq2));
KE = KE1 + KE2 + KE4;
D = simplify(jacobian(jacobian(KE,dq).',dq));

PE = g*m1*p1(2) + g*(m2+m3)*p2(2) + g*m4*p4(2);
G = simplify(jacobian(PE,q).');

syms C
n = length(q);
for k = 1:n
    for j = 1:n
        C(k,j) = 0;
        for i = 1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*dq(i);
        end
    end
end
C = simplify(C);

tau = D*ddq + C*dq + G;

%%
Wheel_R = 0.15; Wheel_M = 0.5;
Bar1_L = 0.3; Bar1_M = 0.1;
Bar2_L = 0.3; Bar2_M = 0.1;
Pen_R = 0.15; Pen_M = 0.5;

[a1,a2,a4] = solve(tau(1)==u1,tau(2)==-u1,tau(3)==u3,ddq1,ddq2,ddq4);

x = [q1 dq1 q2 dq2];
f = [dq1 a1 dq2 a2];
u = [u1 u3];

df_x = jacobian(f,x);
A_h = simplify(subs(df_x,{q1,q2,q3,q4,dq1,dq2,dq3,dq4,u1,u2,u3},{0,0,0,0,...
                                                                 0,0,0,0,...
                                                                 0,0,0}));

df_u = jacobian(f,u);
B_h = simplify(subs(df_u,{q1,q2,q3,q4,dq1,dq2,dq3,dq4,u1,u2,u3},{0,0,0,0,...
                                                                 0,0,0,0,...
                                                                 0,0,0}));

A = eval(subs(A_h,{m1,m2,m3,m4,r1,r2,L1,L2,J1,J2,g},...
                  {Wheel_M*2, Bar1_M, Bar2_M, Pen_M,...
                   Wheel_R, Pen_R, Bar1_L, Bar2_L,...
                   Wheel_M*Wheel_R^2/2, Pen_M*Pen_R^2/2, 9.8}))

B = eval(subs(B_h,{m1,m2,m3,m4,r1,r2,L1,L2,J1,J2,g},...
                  {Wheel_M*2, Bar1_M, Bar2_M, Pen_M,...
                   Wheel_R, Pen_R, Bar1_L, Bar2_L,...
                   Wheel_M*Wheel_R^2/2, Pen_M*Pen_R^2/2, 9.8}))

rank(ctrb(A,B))

% T = 0.02;
% Ad = (eye(length(x))+A*T);
% Bd = B*T;

sdesired = 3*[-1,-2,-3 ,-4];
K1 = place(A,B,sdesired)

Q = [10 0 0 0;
     0 0 0 0;
     0 0 10 0
     0 0 0 1];
R = [5 0;
     0 10];

K2 = lqr(A, B, Q, R)

eig1 = eig(A-B*K1)
eig2 = eig(A-B*K2)
