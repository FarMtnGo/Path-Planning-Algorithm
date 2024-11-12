clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能：五次路径规划，参数显示求解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界条件
qs = 0;   % 起始点
qe = 10;  % 终止点
vs = 0;   % 起始速度
ve = 0;   % 终止速度
as = 0;
ae = 0;
ts = 0;
te = 5;
% 计算规划参数a0,a1,a2,a3
para = QuinticPara(qs, qe, vs, ve, as, ae, ts, te);
i = 1; 
for t = ts: 0.001: te
   time(i) = 0.001*i;
   q(i) = Quintic_position(t, ts, para(1), para(2), para(3), para(4), para(5), para(6));
   qd(i) = Quintic_velocity(t, ts, para(1), para(2), para(3), para(4), para(5), para(6));
   qdd(i) = Quintic_acceleration(t, ts, para(1), para(2), para(3), para(4), para(5), para(6));
   qddd(i) = Quintic_jerk(t, ts, para(1), para(2), para(3), para(4), para(5), para(6));
   i = i + 1;
end
figure(1)
subplot(4, 1, 1)
plot(time, q, 'r', 'LineWidth', 1.5)
grid on
subplot(4, 1, 2)
plot(time, qd, 'b', 'LineWidth', 1.5)
grid on
subplot(4, 1, 3)
plot(time, qdd, 'g', 'LineWidth', 1.5)
grid on
subplot(4, 1, 4)
plot(time, qddd, 'LineWidth', 1.5)
grid on

function para = QuinticPara(qs, qe, vs, ve, as, ae, ts, te)
h = qe-qs;
T = te-ts;
a0 = qs;
a1 = vs;
a2 = as/2;
a3 = (20*h-(8*ve+12*vs)-(3*as-ae)*T^2)/(2*T^3);
a4 = (-30*h-(14*ve+16*vs)*T-(3*as-2*ae)*T^2)/(2*T^4);
a5 = (12*h-(6*ve+6*vs)*T-(as-ae)*T^2)/(2*T^5);
para = [a0,a1,a2,a3,a4,a5];
end

%% 计算位移
function q = Quintic_position(t, ts, a0, a1, a2, a3, a4, a5)
    q = a0 + a1*(t-ts) + a2*(t-ts)^2 + a3*(t-ts)^3 + a4*(t-ts)^4 + a5*(t-ts)^5;
end


%% 计算速度
function qd = Quintic_velocity(t, ts, a0, a1, a2, a3, a4, a5)
    qd = a1 + 2*a2*(t-ts) + 3*a3*(t-ts)^2 + 4*a4*(t-ts)^3 + 5*a5*(t-ts)^4;
end


%% 计算加速度
function qdd = Quintic_acceleration(t, ts, a0, a1, a2, a3, a4, a5)
    qdd = 2*a2 + 6*a3*(t-ts) + 12*a4*(t-ts)^2 + 20*a5*(t-ts)^3;
end

%% 计算加加速度
function qdd = Quintic_jerk(t, ts, a0, a1, a2, a3, a4, a5)
    qdd = 6*a3 + 24*a4*(t-ts) + 60*a5*(t-ts)^2;
end