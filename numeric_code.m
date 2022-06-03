%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 05/2022      %
%  Electricity Numeric Ex   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Contstants: %%%%%%%%

E = 1;              %Electric field const.
B = -1;             %Magnetic field const.
q = 1;              %particle charge const.
m = 1;              %particle mass const.
w = abs(q*B/m);     %Max sequence const.
T = 2*pi/w;         %Max time const.

%% Question 2:

%Set initial values and consts:
t = linspace(0,T);
[~,y,z] = AnaliticalPos(E,B,w,t);
[~,y_speed,z_speed] = AnaliticalSpeed(E,B,w,t);

% PlotFunc(y,z);
xlabel('Y');
ylabel('Z');

% PlotFunc(y_speed,z_speed);
xlabel('V(y)');
ylabel('V(z)');

%% Question 3:
% Use taylor theorem to get r,v on the first-order:
% r = r +v*delta_t.
% v = v +a*delta_t.

TimeJumpsNum = 1000;
delta_t = T/TimeJumpsNum;

%Set values to zero 
y = zeros(TimeJumpsNum,1);
z = zeros(TimeJumpsNum,1);
y_speed = zeros(TimeJumpsNum,1);
z_speed = zeros(TimeJumpsNum,1);

% Get initial speed from analitical answer
[~,y_speed(1),z_speed(1)] = AnaliticalSpeed(E,B,w,0);

%Repeat the taylor sum for every time jump.
for i = 1:TimeJumpsNum
    [pos,speed] = TaylorSum(E,B,m,q,delta_t,[0,y(i),z(i)],[0,y_speed(i),z_speed(i)]);
    y(i+1) = pos(2);
    z(i+1) = pos(3);
    y_speed(i+1) = speed(2);
    z_speed(i+1) = speed(3);
end

%Plot first_order r
% PlotFunc(y,z);
xlabel('Y');
ylabel('Z');

%Plot first_order v
% PlotFunc(y_speed,z_speed);
xlabel('Y Speed');
ylabel('Z Speed');

%% Question 4:

[~,final_y,final_z] = AnaliticalPos(E,B,w,T);

%% For loops parameters
t = 100:10000;
delta_ts = T./t;

%Set values to zero
taylor_error = zeros(length(delta_ts),1); 
midpoint_error = zeros(length(delta_ts),1); 
ronga_kota_error = zeros(length(delta_ts),1); 
rk_k_y = zeros(4,1);
rk_k_z = zeros(4,1);
rk_k_y_speed = zeros(4,1);
rk_k_z_speed = zeros(4,1);
mp_k_y = zeros(2,1);
mp_k_z = zeros(2,1);
mp_k_y_speed = zeros(2,1);
mp_k_z_speed = zeros(2,1);

%% Repeat for different time jumps.
for t = 1:length(delta_ts)
    
    delta_t = delta_ts(t);
    %Set values to zero 
    t_y = zeros(t,1);
    t_z = zeros(t,1);
    t_y_speed = zeros(t,1);
    t_z_speed = zeros(t,1);
    mp_y = zeros(t,1);
    mp_z = zeros(t,1);
    mp_y_speed = zeros(t,1);
    mp_z_speed = zeros(t,1);
    rk_y = zeros(t,1);
    rk_z = zeros(t,1);
    rk_y_speed = zeros(t,1);
    rk_z_speed = zeros(t,1);
    
    % Get initial speed from analitical answer
    [~,t_y_speed(1),t_z_speed(1)] = AnaliticalSpeed(E,B,w,0);
    [~,mp_y_speed(1),mp_z_speed(1)] = AnaliticalSpeed(E,B,w,0);
    [~,rk_y_speed(1),rk_z_speed(1)] = AnaliticalSpeed(E,B,w,0);

    %Repeat the taylor sum for every time jump.
    for i = 1:t

        %%%%%%%% Taylor: %%%%%%%%%
        [pos,speed] = TaylorSum(E,B,m,q,delta_t,[0,t_y(i),t_z(i)],[0,t_y_speed(i),t_z_speed(i)]);
        t_y(i+1) = pos(2);
        t_z(i+1) = pos(3);
        t_y_speed(i+1) = speed(2);
        t_z_speed(i+1) = speed(3);
        
        %%%%%%%% Midpoint: %%%%%%%%
        speed = [0,mp_y_speed(i),mp_z_speed(i)];
        for j = 1:2
            [~,y_acc,z_acc] = AnaliticalAcc(E,B,m,q,[0, speed(2), speed(3)]);
            mp_k_y(j) = delta_t*speed(2);
            mp_k_z(j) = delta_t*speed(3);
            mp_k_y_speed(j) = delta_t*y_acc;
            mp_k_z_speed(j) = delta_t*z_acc;

            [pos,speed] = TaylorSum(E,B,m,q,delta_t/2,[0,mp_k_y(j)/2 + mp_y(i),mp_k_z(j)/2 + mp_z(i)], ...
                [0,mp_k_y_speed(j)/2 + speed(2),mp_k_z_speed(j)/2 + speed(3)]);
        end
        
        mp_y(i+1) = mp_y(i) + mp_k_y(2);
        mp_z(i+1) = mp_z(i) + mp_k_z(2);
        mp_y_speed(i+1) = mp_y_speed(i) + mp_k_y_speed(2);
        mp_z_speed(i+1) = mp_z_speed(i) + mp_k_z_speed(2);        

        %%%%%%%% Ronga-Kota: %%%%%%%%
        speed = [0,rk_y_speed(i),rk_z_speed(i)];
        for j = 1:4
            [~,y_acc,z_acc] = AnaliticalAcc(E,B,m,q,[0, speed(2), speed(3)]);
            rk_k_y(j) = delta_t*speed(2);
            rk_k_z(j) = delta_t*speed(3);
            rk_k_y_speed(j) = delta_t*y_acc;
            rk_k_z_speed(j) = delta_t*z_acc;

            [pos,speed] = TaylorSum(E,B,m,q,delta_t/2,[0,rk_k_y(j)/2 + rk_y(i),rk_k_z(j)/2 + rk_z(i)], ...
                [0,rk_k_y_speed(j)/2 + speed(2),rk_k_z_speed(j)/2 + speed(3)]);
        end
        
        rk_y(i+1) = rk_y(i) + (1/6)*(rk_k_y(1) + 2*rk_k_y(2) + 2*rk_k_y(3) + rk_k_y(4));
        rk_z(i+1) = rk_z(i) + (1/6)*(rk_k_z(1) + 2*rk_k_z(2) + 2*rk_k_z(3) + rk_k_z(4));
        rk_y_speed(i+1) = rk_y_speed(i) + (1/6)*(rk_k_y_speed(1) + 2*rk_k_y_speed(2) + 2*rk_k_y_speed(3) + rk_k_y_speed(4));
        rk_z_speed(i+1) = rk_z_speed(i) + (1/6)*(rk_k_z_speed(1) + 2*rk_k_z_speed(2) + 2*rk_k_z_speed(3) + rk_k_z_speed(4));

    end

    taylor_error(t) = sqrt(abs(final_y - t_y(t+1))^2 + abs(final_z - t_z(t+1))^2);
    midpoint_error(t) = sqrt(abs(final_y - mp_y(t+1))^2 + abs(final_z - mp_z(t+1))^2);
    ronga_kota_error(t) = sqrt(abs(final_y - rk_y(t+1))^2 + abs(final_z - rk_z(t+1))^2);

end

%Plot first order
PlotFunc(delta_ts,taylor_error);
xlabel('delta t');
ylabel('taylor error');
set(gca,'YScale','log');
set(gca,'XScale','log');

%Plot second order
PlotFunc(delta_ts,midpoint_error);
xlabel('delta t');
ylabel('midpoint error');
set(gca,'YScale','log');
set(gca,'XScale','log');

%Plot third order
PlotFunc(delta_ts,ronga_kota_error);
xlabel('delta t');
ylabel('ronga kota error');
set(gca,'YScale','log');
set(gca,'XScale','log');

%% Question 1:
% The particle's movement vs time goes as followes (explained in doc):

function [x,y,z] = AnaliticalPos(e,b,w,t)
    x = 0;
    y = (2*e/(w*b))*cos(w.*t) - 2*e/(w*b);
    z = (2*e/(w*b))*sin(w.*t) + e.*t./b;
end

function [x,y,z] = AnaliticalSpeed(e,b,w,t)
    x = 0;
    y = -2*(e/b)*sin(w*t);
    z = 2*(e/b)*cos(w*t) + e/b;
end

function [x,y,z] = AnaliticalAcc(e,b,m,q,speed_vec)
    %speed_vec = (x_speed,y_speed,z_speed)
    x = 0;
    y = (q/m)*(e-b*speed_vec(3));
    z = (q/m)*b*speed_vec(2);
end

%% Question 3 func:

function [pos,speed] = TaylorSum(e,b,m,q,delta_t,init_pos_vec,init_speed_vec)

    % r = r +v*delta_t.
    pos = init_pos_vec + init_speed_vec.*delta_t;

    % v = v +a*DeltaT.
    [x_acc,y_acc,z_acc] = AnaliticalAcc(e,b,m,q,init_speed_vec);
    speed = init_speed_vec + [x_acc,y_acc,z_acc].*delta_t;
end

%% Misc:
 
%%%%%%%% Functions: %%%%%%%%

function PlotFunc(x,y)
    figure
    hold on
    box on
    grid
    plot(x,y)
    hold off
end 

