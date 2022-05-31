    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 05/2022      %
%  Electricity Numeric Ex   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Contstants: %%%%%%%%

E = 1;              %Electric field const.
B = -1;             %Magnetic field const.
q = 1;              %particle charge const.
m = 1;              %particle mass const.
w = q*B/m;          %Max sequence const.
T = 2*pi/w;         %Max time const.

%% Question 2:

%Set initial values and consts:
t = linspace(0,T);
[~,y,z] = AnaliticalPos(E,B,w,t);
[~,y_speed,z_speed] = AnaliticalSpeed(E,B,w,t);

PlotFunc(y,z);
xlabel('Y');
ylabel('Z');

PlotFunc(y_speed,z_speed);
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

[~,y_speed(1),z_speed(1)] = AnaliticalSpeed(E,B,w,0);

for i = 1:TimeJumpsNum
    [pos,speed] = TaylorSum(E,B,m,q,delta_t,[0,y(i),z(i)],[0,y_speed(i),z_speed(i)]);
    y(i+1) = pos(2);
    z(i+1) = pos(3);
    y_speed(i+1) = speed(2);
    z_speed(i+1) = speed(3);
    
end

%Plot first_order r
PlotFunc(y,z);
xlabel('Y');
ylabel('Z');

%Plot first_order v
PlotFunc(y_speed,z_speed);
xlabel('Y Speed');
ylabel('Z Speed');

%% Question 4:

%%%%%%%% midpoint: %%%%%%%%

% NumOfTimeJumps = 1000;

for i = 1:TimeJumpsNum

    k1 = delta_t*y_speed(i);

    % r = r +v*delta_t.
    y(i+1) = y(i) + y_speed(i)*delta_t;
    z(i+1) = z(i) + z_speed(i)*delta_t;

    % v = v +a*DeltaT.
    [~, y_acc(i),z_acc(i)] = AnaliticalAcc(E,B,m,q,[0,y_speed(i),z_speed(i)]);
    y_speed(i+1) = y_speed(i) + y_acc(i)*delta_t;
    z_speed(i+1) = z_speed(i) + z_acc(i)*delta_t;
    
end

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

%% Question 2 func:

function [pos,speed] = TaylorSum(e,b,m,q,delta_t,init_pos_vec,init_speed_vec)

    % r = r +v*delta_t.
    pos = init_pos_vec + init_speed_vec.*delta_t;

    % v = v +a*DeltaT.
    [x_acc,y_acc,z_acc] = AnaliticalAcc(e,b,m,q,init_speed_vec);
    speed = init_speed_vec + [x_acc,y_acc,z_acc].*delta_t;
end

%% Misc:
 
%%%%%%%% Functions: %%%%%%%%

function p = PlotFunc(x,y)
    figure
    hold on
    box on
    grid
    plot(x,y)
    hold off
end 

