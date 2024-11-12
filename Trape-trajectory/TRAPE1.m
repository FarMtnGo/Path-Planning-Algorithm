clc;clear;close all;
%% 梯形路径规划
ps = 0;     % 初始值
pe = 15;   % 终点值
Vel = 2*0.2;    % 最大速度
Acc = 0.5*0.2;  % 加速度
Dec = 0.5*0.2;  % 减速度
dt = 0.02;  % 规划时间间隔

[t1,t2,t3] = Get_Trape_Para(Vel,Acc,Dec,abs(pe-ps));
ts = t1 + t2 + t3;

[qd,dqd,ddqd] = Get_Trape_Plan(t1,t2,t3,Acc,Dec,ps,sign(pe-ps),dt,ts);
tsm = linspace(0,ts,ceil(ts/dt));

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

function [P_t1,P_t2,P_t3] = Get_Trape_Para(P_Vel,P_Acc,P_Dec,sq)
    ta = P_Vel/P_Acc;
    td = P_Vel/P_Dec;
    sa = 0.5*P_Acc*ta^2;
    sd = 0.5*P_Dec*td^2;
    if sa+sd > sq
        P_Vel = sqrt(2*sq/(1/P_Acc+1/P_Dec));
        P_t1 = P_Vel/P_Acc;
        P_t2 = 0;
        P_t3 = P_Vel/P_Dec;
    else
        P_t1 = ta;
        P_t2 = (sq-sa-sd)/P_Vel;
        P_t3 = td;
    end
end

function [qd,dqd,ddqd] = Get_Trape_Plan(P_t1,P_t2,P_t3,P_Acc,P_Dec,p_s,A_Sign,dt,ts)
    data_num = ceil(ts/dt);
    P_Acc = A_Sign*P_Acc;
    P_Dec = A_Sign*P_Dec;
    V1 = P_t1*P_Acc;
    V2 = V1;
    P1 = p_s + 0.5*P_Acc*P_t1^2;
    P2 = P1 + V1*P_t2;
    for i = 1:data_num
        if (i*dt-P_t1) <= 0.00001
            temp_t = i*dt;
            ddqd(i) = P_Acc;
            dqd(i) = P_Acc*temp_t;
            qd(i) = p_s + 0.5*P_Acc*temp_t^2;
        elseif (i*dt - P_t1 - P_t2) <= 0.00001
            temp_t = i*dt - P_t1;
            ddqd(i) = 0;
            dqd(i) = V1;
            qd(i) = P1 + V1*temp_t;
        elseif (i*dt - P_t1 - P_t2 - P_t3) <= 0.00001
            temp_t = i*dt - P_t1 - P_t2;
            ddqd(i) = -P_Dec;
            dqd(i) = V2 - P_Dec*temp_t;
            qd(i) = P2 + 0.5*(V2 + dqd(i))*temp_t;
        else
            temp_t = P_t3;
            ddqd(i) = 0;
            dqd(i) = 0;
            qd(i) = P2 + 0.5*(V2 + dqd(i))*temp_t;
        end
    end
end