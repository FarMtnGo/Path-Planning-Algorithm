clc;clear;close all;
% 边界条件
qs = 10;   % 起始点
qe = 0;    % 终止点
vs = 0;    % 起始速度
ve = 0;    % 终止速度
ts = 0;
te = 5;
% 计算规划参数a0,a1,a2,a3
para = CubicPara(qs, qe, vs, ve, ts, te);
i = 1; 
for t = ts: 0.001: te
   time(i) = 0.001*i;
   q(i) = Cubic_position(t, ts, para(1), para(2), para(3), para(4));
   qd(i) = Cubic_velocity(t, ts, para(1), para(2), para(3), para(4));
   qdd(i) = Cubic_acceleration(t, ts, para(1), para(2), para(3), para(4));
   i = i + 1;
end
figure(1)
subplot(3, 1, 1)
plot(time, q, 'r', 'LineWidth', 1.5)
grid on
xlabel('t(s)'); ylabel('pos');
subplot(3, 1, 2)
plot(time, qd, 'b', 'LineWidth', 1.5)
grid on
xlabel('t(s)'); ylabel('vel');
subplot(3, 1, 3)
plot(time, qdd, 'g', 'LineWidth', 1.5)
grid on
xlabel('t(s)'); ylabel('acc');

function para = CubicPara(qs, qe, vs, ve, ts, te)
h = qe-qs;
T = te-ts;
a0 = qs;
a1 = vs;
a2 = 3*h/T^2 - (2*vs+ve)/T;
a3 = -2*h/T^3 + (vs+ve)/T^2;
para = [a0,a1,a2,a3];
end

%% 计算位移
function q = Cubic_position(t, ts, a0, a1, a2, a3)
    q = a0 + a1*(t-ts) + a2*(t-ts)^2 + a3*(t-ts)^3;
end


%% 计算速度
function qd = Cubic_velocity(t, ts, a0, a1, a2, a3)
    qd = a1 + 2*a2*(t-ts) + 3*a3*(t-ts)^2;
end


%% 计算加速度
function qdd = Cubic_acceleration(t, ts, a0, a1, a2, a3)
    qdd = 2*a2 + 6*a3*(t-ts);
end
