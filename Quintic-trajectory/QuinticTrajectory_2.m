clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能：五次路径规划，并考虑对中间点的速度和加速度
%       进行连续化处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 边界条件
ts=[0,2,4,8,10];
qs=[10,20,0,30,40];
vs = 0;
ve = 0;
as = 0;
ae = 0;
%起点和终点加速度假设为0，中间点加速度都初始化为0，

dt = 0.001;
nums = length(ts);
% 处理中间速度
for i=1:nums
    if i==1
        Vec_vs(i) = vs;    % 起始
    elseif i==nums
        Vec_vs(i) = ve;    % 终止
    else
        dv1 = (qs(i)-qs(i-1))/(ts(i)-ts(i-1));  % 中间
        dv2 = (qs(i+1)-qs(i))/(ts(i+1)-ts(i));
        if ( dv1 >=0 && dv2 >=0 ) || ( dv1 <=0 && dv2 <=0 )
            Vec_vs(i) = (dv1 + dv2)/2;          % 速度同向取均值
        else
            Vec_vs(i) = 0;                      % 速度反向取0值
        end
    end
end
% 处理中间加速度
for i=1:nums
    if i==1
        Vec_as(i) = as;
    elseif i==nums
        Vec_as(i) = ae;
    else
        dv1 = (Vec_vs(i)-Vec_vs(i-1))/(ts(i)-ts(i-1));
        dv2 = (Vec_vs(i+1)-Vec_vs(i))/(ts(i+1)-ts(i));
        if ( dv1 >=0 && dv2 >=0 ) || ( dv1 <=0 && dv2 <=0 )
            Vec_as(i) = (dv1 + dv2)/2;
        else
            Vec_as(i) = 0;
        end
    end
end
j = 1;
ti = ts(1);
% 开始规划
while ti <= ts(end)
    for i=1:nums-1
        if ti >= ts(i) && ti < ts(i+1)
            para = Get_quintic_Para(qs(i), qs(i+1), Vec_vs(i), Vec_vs(i+1), Vec_as(i), Vec_as(i+1), ts(i), ts(i+1));
            [qi,dqi,ddqi] = Get_quintic_Plan(ti,para);
        end
    end
    qd(j) = qi;
    dqd(j) = dqi;
    ddqd(j) = ddqi;
    tsm(j) = ti;
    j = j + 1;
    ti = ti + dt;
end

figure(1);
subplot(3,1,1)
plot(tsm,qd,'Color','r','LineStyle','-','LineWidth',2);
grid on;
xlabel('Time(s)')
ylabel('q')
subplot(3,1,2)
plot(tsm,dqd,'Color','g','LineStyle','-','LineWidth',2);
grid on;
xlabel('Time(s)')
ylabel('dq')
subplot(3,1,3)
plot(tsm,ddqd,'Color','b','LineStyle','-','LineWidth',2);
grid on;
xlabel('Time(s)')
ylabel('ddq')

function para = Get_quintic_Para(q0,qf,v0,vf,a0,af,t0,tf)
    A = zeros(6,6);
    A(1,1) = 1; A(1,2) = t0; A(1,3) = t0^2; A(1,4) = t0^3; A(1,5) = t0^4; A(1,6) = t0^5;
    A(2,1) = 1; A(2,2) = tf; A(2,3) = tf^2; A(2,4) = tf^3; A(2,5) = tf^4; A(2,6) = tf^5;
    A(3,1) = 0; A(3,2) = 1; A(3,3) = 2*t0; A(3,4) = 3*t0^2; A(3,5) = 4*t0^3; A(3,6) = 5*t0^4;
    A(4,1) = 0; A(4,2) = 1; A(4,3) = 2*tf; A(4,4) = 3*tf^2; A(4,5) = 4*tf^3; A(4,6) = 5*tf^4;
    A(5,1) = 0; A(5,2) = 0; A(5,3) = 2; A(5,4) = 6*t0; A(5,5) = 12*t0^2; A(5,6) = 20*t0^3;
    A(6,1) = 0; A(6,2) = 0; A(6,3) = 2; A(6,4) = 6*tf; A(6,5) = 12*tf^2; A(6,6) = 20*tf^3;
    B = [q0; qf; v0; vf; a0; af];
    para = A\B;
end

function [qi,dqi,ddqi] = Get_quintic_Plan(ti,para)
    a0 = para(1);
    a1 = para(2);
    a2 = para(3);
    a3 = para(4);
    a4 = para(5);
    a5 = para(6);
    qi = a0 + a1*ti + a2*ti^2 + a3*ti^3 + a4*ti^4 + a5*ti^5;
    dqi = a1 + 2*a2*ti + 3*a3*ti^2 + 4*a4*ti^3 + 5*a5*ti^4;
    ddqi = 2*a2 + 6*a3*ti + 12*a4*ti^2 + 20*a5*ti^3;
end