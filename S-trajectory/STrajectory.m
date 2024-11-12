%% S曲线规划
% 边界条件
qs = 10;   % 起始点
qe = 0;  % 终止点
vs = 0;   % 起始速度
ve = 0;   % 终止速度
vmax = 5; % 最大速度
amax = 5; % 最大加速度
jmax = 20; % 最大加加速度
sigma = sign(qe - qs);
% 计算规划参数Ta, Tv, Td, Tja, Tjd, qs, qe, vs, ve, vlim, amax, amin, alima, alimd, jmax, jmin
para = STrajectoryPara(qs, qe, vs, ve, vmax, amax, jmax);
i = 1; 
T = para(1) + para(2) + para(3)
for t = 0: 0.001: T
   time(i) = 0.001*i;
   q(i) = S_position(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8), para(9), para(10), para(11), para(12), para(13), para(14), para(15), para(16));
   qd(i) = S_velocity(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8), para(9), para(10), para(11), para(12), para(13), para(14), para(15), para(16));
   qdd(i) = S_acceleration(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8), para(9), para(10), para(11), para(12), para(13), para(14), para(15), para(16));
   qddd(i) = S_jerk(t, para(1), para(2), para(3), para(4), para(5), para(6), para(7), para(8), para(9), para(10), para(11), para(12), para(13), para(14), para(15), para(16));
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
%% S曲线参数计算（S型速度规划，又称七段式轨迹）
function para = STrajectoryPara(qs, qe, vs, ve, vmax, amax, jmax)
% 得到规划参数Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, amax, amin, alima, alimd, jmax, jmin
% 输入：
% 边界条件qs, qe, vs, ve
% 约束条件vmax, amax, jamx
% 输出：
% 加速时间：Ta
% 匀速时间：Tv
% 减速时间：Td
% 加加速时间：Tja
% 减加速时间：Tjd

% 默认条件：
vmin = -vmax; amin = -amax; jmin = -jmax;

%% 利用公式（3.31）（3.32）转化得到实际的q_0、q_1、v_max、a_max
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

%% 判断是否达到最大速度
if ((v_max - v_s)*j_max < a_max^2) 
    Tja = sqrt((v_max - v_s) / j_max); % 达不到a_max
    Ta = 2*Tja;
    a_lima = j_max * Tja;
else
    Tja = a_max / j_max; % 能够达到a_max
    Ta = Tja + (v_max - v_s) / a_max;
    a_lima = a_max;
end
if ((v_max - v_e)*j_max < a_max^2)
    Tjd = sqrt((v_max - v_e) / j_max); % 达不到a_min
    Td = 2*Tjd;
    a_limd =  -j_max * Tjd;
else
    Tjd = a_max / j_max; % 能够达到a_min
    Td = Tjd + (v_max - v_e) / a_max;
    a_limd = -a_max;
end
% 根据（3.25）计算匀速段时间
Tv = (q_e - q_s)/v_max - (Ta/2)*(1 + v_s/v_max) - (Td/2)*(1 + v_e/v_max);

%% 对Tv进行讨论
if (Tv > 0)
    % 达到最大速度v_max，即存在匀速阶段
    vlim = v_max;
    T = Ta + Tv + Td;
    para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
    return;
else
    % 达不到最大速度，即匀速阶段Tv=0
    % 假设最大加速度和最小加速度均能达到
    Tv = 0;
    Tj = a_max / j_max;
    Tja = Tj;
    Tjd = Tj;
    delta = (a_max^4/j_max^2) + 2*(v_s^2 + v_e^2) + a_max*(4*(q_e - q_s) - 2*(a_max/j_max)*(v_s + v_e));
    Ta = ((power(a_max, 2)/j_max) - 2.0*v_s + sqrt(delta)) / (2.0*a_max);
    Td = ((power(a_max, 2)/j_max) - 2.0*v_e + sqrt(delta)) / (2.0*a_max);
    % 对Ta和Td进行讨论
    if (Ta < 0 || Td < 0)
        if (Ta < 0)
            % 没有加速段，只有减速段
            Ta = 0; Tja = 0;
            Td = 2*(q_e - q_s) / (v_s + v_e);
            Tjd = (j_max*(q_e - q_s) - sqrt(j_max*(j_max*power(q_e - q_s, 2) + power(v_e + v_s, 2)*(v_e - v_s)))) / (j_max*(v_e + v_s));
            a_lima = 0;
            a_limd = -j_max*Tjd;
            vlim = vs;
            para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
            return;
        elseif (Td < 0)
            % 没有减速段，只有加速段
            Td = 0; Tjd = 0;
            Ta = 2*(q_e - q_s) / (v_s + v_e);
            Tja = (j_max*(q_e - q_s) - sqrt(j_max*(j_max*power(q_e - q_s, 2)) - power(v_e + v_s, 2)*(v_e - v_s))) / (j_max*(v_e + v_s));
            a_lima = j_max*Tja;
            a_limd = 0;
            vlim = v_s + a_lima*(Ta - Tja);
            para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
            return;
        end
    elseif (Ta >= 2*Tj && Td >= 2*Tj)
        % 加速段和减速段都能达到最大加速度
        a_lima = a_max;
        a_limd = -a_max;
        vlim = vs + a_lima*(Ta - Tj);
        para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
        return;
    else
        % 加速和减速阶段至少有一段不能达到最大加速度
        lambda = 0.99; % 系统取0<lambda<1
        while (Ta < 2*Tj || Td < 2*Tj)
            % 循环
            a_max = lambda*a_max;
            Tv = 0;
            Tj = a_max / j_max;
            Tja = Tj;
            Tjd = Tj;
            delta = (a_max^4/j_max^2) + 2*(v_s^2 + v_e^2) + a_max*(4*(q_e - q_s) - 2*(a_max/j_max)*(v_s + v_e));
            Ta = ((power(a_max, 2)/j_max) - 2.0*v_s + sqrt(delta)) / (2.0*a_max);
            Td = ((power(a_max, 2)/j_max) - 2.0*v_e + sqrt(delta)) / (2.0*a_max);
            if (Ta < 0 || Td < 0)
                if (Ta < 0)
                    % 没有加速段，只有减速段
                    Ta = 0; Tja = 0;
                    Td = 2*(q_e - q_s) / (v_s + v_e);
                    Tjd = (j_max*(q_e - q_s) - sqrt(j_max*(j_max*power(q_e - q_s, 2) + power(v_e + v_s, 2)*(v_e - v_s)))) / (j_max*(v_e + v_s));
                    a_lima = 0;
                    a_limd = -j_max*Tjd;
                    vlim = vs;
                    para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
                    return;
                elseif (Td < 0)
                    % 没有减速段，只有加速段
                    Td = 0; Tjd = 0;
                    Ta = 2*(q_e - q_s) / (v_s + v_e);
                    Tja = (j_max*(q_e - q_s) - sqrt(j_max*(j_max*power(q_e - q_s, 2)) - power(v_e + v_s, 2)*(v_e - v_s))) / (j_max*(v_e + v_s));
                    a_lima = j_max*Tja;
                    a_limd = 0;
                    vlim = v_s + a_lima*(Ta - Tja);
                    para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
                    return;
                end
            elseif (Ta >= 2*Tj && Td >= 2*Tj)
                % 加速段和减速段都能达到最大加速度
                a_lima = a_max;
                a_limd = -a_max;
                vlim = vs + a_lima*(Ta - Tj);
                para = [Ta, Tv, Td, Tja, Tjd, q_s, q_e, v_s, v_e, vlim, a_max, a_min, a_lima, a_limd, j_max, j_min];
                return;
            end
        end
    end
end
end

%% 计算位移
function q = S_position(t, Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, amax, amin, alima, alimd, jmax, jmin)
T = Ta + Tv + Td;
% 加速段
if (t >= 0 && t < Tj1)
    q = q0 + v0*t + jmax*t^3/6;
elseif (t >= Tj1 && t < Ta - Tj1)
    q = q0 + v0*t +(alima/6)*(3*t^2 - 3*Tj1*t + Tj1^2);
elseif (t >= Ta - Tj1 && t < Ta)
    q = q0 + (vlim + v0)*(Ta/2) - vlim*(Ta - t) - jmin*((Ta - t)^3/6);
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    q = q0 + (vlim + v0)*(Ta/2) + vlim*(t - Ta);
% 减速段
elseif (t >= Ta + Tv && t < T - Td + Tj2)
    q = q1 - (vlim + v1)*(Td/2) + vlim*(t - T + Td) - jmax*(power(t - T + Td, 3)/6);
elseif (t >= T - Td + Tj2 && t < T - Tj2)
    q = q1 - (vlim + v1)*(Td/2) + vlim*(t - T + Td) + (alimd/6)*(3*power(t - T + Td, 2) - 3*Tj2*(t - T + Td) + Tj2^2);
elseif (t >= T - Tj2 && t <= T)
    q = q1 - v1*(T - t) - jmax*(power(T - t, 3)/6);
end
end


%% 计算速度
function qd = S_velocity(t, Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, amax, amin, alima, alimd, jmax, jmin)
T = Ta + Tv + Td;
if (t >= 0 && t < Tj1)
    qd = v0 + jmax*(t^2/2);
elseif (t >= Tj1 && t < Ta - Tj1)
    qd = v0 + alima*(t - Tj1/2);
elseif (t >= Ta - Tj1 && t < Ta)
    qd = vlim + jmin*(power(Ta - t, 2)/2);
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qd = vlim;
% 减速段
elseif (t >= Ta + Tv && t < T - Td + Tj2)
    qd = vlim - jmax*(power(t - T + Td, 2)/2);
elseif (t >= T - Td + Tj2 && t < T - Tj2)
    qd = vlim + alimd*(t - T + Td - Tj2/2);
elseif (t >= T - Tj2 && t <= T)
    qd = v1 + jmax*(power(t - T, 2)/2);
end
end


%% 计算加速度
function qdd = S_acceleration(t, Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, amax, amin, alima, alimd, jmax, jmin)
T = Ta + Tv + Td;
if (t >= 0 && t < Tj1)
    qdd = jmax*t;
elseif (t >= Tj1 && t < Ta - Tj1)
    qdd = alima;
elseif (t >= Ta - Tj1 && t < Ta)
    qdd = -jmin*(Ta - t);
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qdd = 0;
% 减速段
elseif (t >= Ta + Tv && t < T - Td + Tj2)
    qdd = -jmax*(t - T + Td);
elseif (t >= T - Td + Tj2 && t < T - Tj2)
    qdd = alimd;
elseif (t >= T - Tj2 && t <= T)
    qdd = -jmax*(T - t);
end
end


%% 计算加加速度
function qddd = S_jerk(t, Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, amax, amin, alima, alimd, jmax, jmin)
T = Ta + Tv + Td;
if (t >= 0 && t < Tj1)
    qddd = jmax;
elseif (t >= Tj1 && t < Ta - Tj1)
    qddd = 0;
elseif (t >= Ta - Tj1 && t < Ta)
    qddd = jmin;
% 匀速段
elseif (t >= Ta && t < Ta + Tv)
    qddd = 0;
% 减速段
elseif (t >= Ta + Tv && t < T - Td + Tj2)
    qddd = -jmax;
elseif (t >= T - Td + Tj2 && t < T - Tj2)
    qddd = 0;
elseif (t >= T - Tj2 && t <= T)
    qddd = jmax;
end
end
