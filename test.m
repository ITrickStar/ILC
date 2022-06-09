addpath(genpath('D:\packages\sedumi'));
addpath(genpath('D:\packages\YALMIP'));

Ts = 0.01;
Time = 1;
T = Time/Ts;
d = 3;
N = 30;

% параметры передаточной функции
y =@(s) 23.7356*(s+661.2);
x =@(s) s*(s^2+426.7*s+1.744*10^5);

X = zeros(T+1,1);
for t=d+1:T+1
    X(t,1) = x(t-d);
end

Y = zeros(T+1,1);
for t=1:T+1
    Y(t,1) = y(t);
end

% дискретная модель
sys = ss(x, y);
sysd = c2d(sys,Ts);
A = sysd.a;
B = sysd.b;
C = sysd.c;

I = eye(1);
A_=[ A; 0; -CA; I];
B_ = [B; -C*B];
C_ = [C, 0; 0, I];

Q = diag(1, 1, 1, 50);
R = 10-3;
K = [-962.2, 185.27];

% матричное неравенство
Cond11 = X;
Cond12 = (X*A_+B_*Y*C_)';
Cond13 = X;
Cond14 = (Y*C_)';
Cond21 = A_*X+B_*Y*C_;
Cond22 = X;
Cond31 = X;
Cond33 = inv(Q);
Cond41 = YC_;
Cond44 = inv(R);
Cond = [Cond11, Cond12, Cond13, Cond14; ...
        Cond21, Cond22, 0,      0; ...
        Cond31, 0,      Cond33, 0 ...
        Cond41, 0,      0,      Cond44];
options = sdpsettings('solver', 'sedumi', 'verbose');
solved = solvesdp(Cond >= 0, options);

% начальные значения
x0 = [0;0;0;0];
xhat0 = [0;0;0;0];

for k=1:N+1
    x{1,k}=x0;
    xhat{1,k}=xhat0;
    y(1,t)=C*x0;
end

for t=1:T+1
    u(t,1)=0;
    x{t+1,1}=A*x{t,1};
    xhat{t+1,1}=A*xhat{t,1};
    y(t,1)=C*x{t,1};
end

for k=2:N+1
    for t=1:T
        u(t,k)=u(t,k-1)+K(1)*(xhat{t,k}-xhat{t,k-1})+K(2)*(yref(t+1,1)-C*x{t+1,k-1});

        x{t+1,k}=A*x{t,k}+B*u(t,k);
        y(t,k)=C*x{t,k};

        xhat{t+1,k}=A*xhat{t,k}+B*u(t,k)-C*xhat{t,k};
        e(t,k-1)= Y(t,1)-y(t,k-1);
        E(1,k-1) = (E(1,k-1) + (e(t,k-1)^2));
    end
    E(1,k-1) = sqrt(E(1,k-1)/T); % среднеквадратическая ошибка обучения
end

% график
[XX, YY] = meshgrid(0:T-1,0:N-1);
figure(1);
plot(YY, E', 'b','LineWidth', 2);
title(str, 'FontSize',14);
xlabel('k','FontSize',14);
ylabel('E','FontSize',14);
grid on;
