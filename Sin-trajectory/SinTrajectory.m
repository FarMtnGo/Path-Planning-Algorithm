%% 正弦平方速度规划
% 边界条件
qs = 10;   % 起始点
qe = 0;  % 终止点
vs = 0;   % 起始速度
ve = 0;   % 终止速度
vmax = 5; % 最大速度
amax = 5; % 最大加速度
jmax = 20; % 最大加加速度
sigma = sign(qe - qs);
% 计算规划参数Ta, Tv, Td, qs, qe, vs, ve, vlim, amax, amin, jmax, jmin
para = SinTrajectoryPara(qs, qe, vs, ve, vmax, amax, jmax);
i = 1; 
T = para(1) + para(2) + para(3)
for t = 0: 0.001: T
   time(i) = 0.001*i;
   q(i) = Sin_position(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8));
   qd(i) = Sin_velocity(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8));
   qdd(i) = Sin_acceleration(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8));
   qddd(i) = Sin_jerk(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8));
   i = i + 1;
end
q = sigma*q;
qd = sigma*qd;
qdd = sigma*qdd;
qddd = sigma*qddd;
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

%% 正弦曲线参数计算
function para = SinTrajectoryPara(qs, qe, vs, ve, vmax, amax, jmax)
% 得到规划参数Ta, Tv, Td, q0, q1, v0, v1, vlim, amax, amin, jmax, jmin
% 输入：
% 边界条件qs, qe, vs, ve
% 约束条件vmax, amax, jamx
% 输出：
% 加速时间：Ta
% 匀速时间：Tv
% 减速时间：Td

% 默认条件：
vmin = -vmax; amin = -amax; jmin = -jmax;
%% 实际的q_0、q_1、v_max、a_max
sigma = sign(qe - qs);  
q_s = sigma*qs;
q_e = sigma*qe;
v_s = sigma*vs;
v_e = sigma*ve;
v_max = ((sigma+1)/2)*vmax + ((sigma-1)/2)*vmin;
v_min = ((sigma+1)/2)*vmin + ((sigma-1)/2)*vmax;
a_max = ((sigma+1)/2)*amax + ((sigma-1)/2)*amin;
a_min = ((sigma+1)/2)*amin + ((sigma-1)/2)*amax;
j_max = ((sigma+1)/2)*jmax + ((sigma-1)/2)*jmin;
j_min = ((sigma+1)/2)*jmin + ((sigma-1)/2)*jmax;
% 加速段时间
Ta = max(2*(v_max-v_s)/a_max, sqrt(2*pi*(v_max-v_s)/j_max));
% 加速段距离
La = (v_s+v_max)*Ta/2;
% 减速段时间
Td = max(2*(v_max-v_e)/a_max, sqrt(2*pi*(v_max-v_e)/j_max));
% 减速段距离
Ld = (v_e+v_max)*Td/2;
% 匀速段距离
Lv = (q_e - q_s)-La-Ld;
Tv = Lv/v_max;
if Tv >= 0
    % 上面假设成立，完成参数计算
    vlim = v_max;
    para = [Ta, Tv, Td, q_s, q_e, v_s, v_e, vlim, a_max, a_min, j_max, j_min];
    return;
else
    % 上面假设不成立，需要重新计算实际规划最大速度vlim
    Tv = 0;
    fun = @(x) SinVlim(x,q_s,q_e,v_s,v_e,a_max,j_max);
    vlim = fzero(fun,v_max);
    Ta = max(2*(vlim-v_s)/a_max, sqrt(2*pi*(vlim-v_s)/j_max));
    La = (v_s+vlim)*Ta/2;
    
    Td = max(2*(vlim-v_e)/a_max, sqrt(2*pi*(vlim-v_e)/j_max));
    Ld = (v_e+vlim)*Td/2;
    para = [Ta, Tv, Td, q_s, q_e, v_s, v_e, vlim, a_max, a_min, j_max, j_min];
    return;
end
end

function y = SinVlim(x,qs,qe,vs,ve,amax,jmax)
Ta = max(2*(x-vs)/amax, sqrt(2*pi*(x-vs)/jmax));
La = (vs+x)*Ta/2;

Td = max(2*(x-ve)/amax, sqrt(2*pi*(x-ve)/jmax));
Ld = (ve+x)*Td/2;

y = (qe - qs)-La-Ld;
end

%% 计算位移
function q = Sin_position(t, Ta, Tv, Td, q0, q1, v0, v1, vlim)
T = Ta + Tv + Td;
% 加速段
if (t >= 0 && t < Ta)
    q = q0 + v0*t + (vlim-v0)*(Ta*(cos(2*pi*t/Ta)-1)/(2*pi)+pi*t^2/Ta)/(2*pi);
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    q = (v0+vlim)*Ta/2 + vlim*(t - Ta);
% 减速段
elseif (t >= Ta + Tv && t <= T)
    q = q1 - v1*(T - t) - (vlim-v1)*(Td*(cos(2*pi*(T-t)/Td)-1)/(2*pi)+pi*(T-t)^2/Td)/(2*pi);

end
end


%% 计算速度
function qd = Sin_velocity(t, Ta, Tv, Td, q0, q1, v0, v1, vlim)
T = Ta + Tv + Td;
if (t >= 0 && t < Ta)
    qd = v0 - (vlim-v0)*(sin(2*pi*t/Ta)-2*pi*t/Ta)/(2*pi);
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qd = vlim;
% 减速段
elseif (t >= Ta + Tv && t <= T)
    qd = v1 - (vlim-v1)*(sin(2*pi*(T-t)/Td)-2*pi*(T-t)/Td)/(2*pi);
end
end


%% 计算加速度
function qdd = Sin_acceleration(t, Ta, Tv, Td, q0, q1, v0, v1, vlim)
T = Ta + Tv + Td;
if (t >= 0 && t < Ta)
    qdd = 2*(vlim-v0)*sin(pi*t/Ta)^2/Ta;
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qdd = 0;
% 减速段
elseif (t >= Ta + Tv && t <= T)
    qdd = -2*(vlim-v1)*sin(pi*(T-t)/Td)^2/Td;
end
end

%% 计算加加速度
function qddd = Sin_jerk(t, Ta, Tv, Td, q0, q1, v0, v1, vlim)
T = Ta + Tv + Td;
if (t >= 0 && t < Ta)
    qddd = 2*pi*(vlim-v0)*sin(2*pi*t/Ta)/Ta^2;
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qddd = 0;
% 减速段
elseif (t >= Ta + Tv && t <= T)
    qddd = 2*pi*(vlim-v1)*sin(2*pi*(T-t)/Td)/Td^2;
end
end