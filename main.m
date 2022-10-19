clear all;clc
z = tf('z');
Gz = tf([0 0.08123],[1 -0.9895],0.01);
[A,B,C,D] = tf2ss([0 0.08123],[1 -0.9895]);
sys = ss(A,B,C,D,0.01);
N = 100;
u = 0;
gamma = 0.01;
%%  Apply Disturbance 
w = wgn(N,1,1);
v = w;
% Rww = 10;
% Rvv = 10;
% w = unifrnd(1,Rww,N,1);
% v = unifrnd(1,Rvv,N,1);

%%  Nominal model
x_nom = .1*ones(1,size(A,1));
for i = 1:N
    x_nom(:,i+1) = A*x_nom(:,i)+B*u;
    y_nom(i) = C*x_nom(:,i);
end
%% Inject Fault 
time = [30 45 65]; 
% x = .1*ones(1,time(1),3);
x = .1*ones(1,N,3);
y = zeros(N,1);
tout = 80;
for j = 1:3
for i = time(j):tout
    x(:,i+1,j) = A*x(:,i,j)+B*u+w(i);
    y(i,j) = C*x(:,i,j)+v(i);
end
end
%% Kalman Conventional 
Rww = .01;
Rvv = 1;
[Kest,L,S] = kalman(sys,N,Rww,Rvv);
% y_hat_kf(1:N) = Cf*x_hat_kf(1:N);
%% Estimated state 
x_hat = 0.1*ones(1,N,3);
for j = 1:3
for i = 1:N
    x_hat(:,i+1,j) = A*x_hat(:,i,j)+L*(y(i,j)-C*x_hat(:,i,j));
    y_hat(i,j) = C*x_hat(:,i,j);
end
end
%% Residual signal 
for i = 1:3
    r(:,i) = (y(:,i)-y_hat(:,i));
end
%% Residual Evaluation
% gamma = 0.1;
count = zeros(1,3);
fault = zeros(1,N);
for j = 1:3
for i = 1:N
    if r(i,j)>gamma||r(i,j)<-gamma
        count(j) = count(j)+1; 
        fault(j,i) = 1;
    end
    
end
% fault(j,N-count(j):N) = 1;
end
%% Fault Tolerant Control using VI 
Q = .01;
R = 10;
[P,K,G] = value_iteration(A,N,B,Q,R);
for j = 1:3
    x_hat_new(1,:,j) = x_hat(1,:,j);
for i=1:N
    u(i,j) = -K(i)*x_hat_new(1,i,j);
    x_hat_new(1,i+1,j) = A*x_hat(1,i,j)+B*u(i,j);
end
end
%%
% fprintf('Q = %f\n',Q)
% fprintf('R = %f\n',R)
% fprintf('Norm K : %f \n',norm(K))
% fprintf('Norm u : %f %f %f\n',norm(u(:,1)),norm(u(:,2)),norm(u(:,3)))
fprintf('Rww = %f\n',Rww)
fprintf('Rvv = %f\n',Rvv)
fprintf('Norm L : %f \n',norm(L))
fprintf('Norm r : %f %f %f\n',norm(r(:,1)),norm(r(:,2)),norm(r(:,3)))

%% Value Function
xx = x_hat_new(1,1:N,1);
J = value_func(xx,u(:,1)',Q,R,N);
% %% Plot
figure(1);clf
for i = 1:3
subplot(3,1,i)
plot(r(:,i),'k','LineWidth',2)
hold on
plot(fault(i,:),'b','LineWidth',2)
grid on
xlabel('Iteration step')
title("Fault injected on " + time(i)+ "th secs until 80th secs",'Interpreter','latex');
hold on
legend('Residual Signals','Residual Evaluation','location','NorthWest')
end
figure(2);clf
for i = 1:3
subplot(3,1,i)
plot(r(:,i),'k','LineWidth',2)
hold on
plot(u(:,i),'m','LineWidth',2)
grid on
xlabel('Iteration step')
title("Fault injected on " + time(i)+ "th secs until 80th secs",'Interpreter','latex');
hold on
legend('Residual Signals','Control Signal','location','NorthWest')
end
