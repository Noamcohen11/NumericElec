%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 05/2022      %
%  Electricity Numeric Ex   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Contstants: %%%%%%%%

E = 1;                  %Electric field const.
B = 1;                  %Magnetic field const.
q = 1;                  %particle charge const.
m = 1;                  %particle mass const.
w = q*B/m;              %Max sequence const.
T = abs(2*pi/w);        %Max time const.

QRun = [1 0 0 0 0 0 1];

%% Question 2:
if QRun(2) 
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
end

%% Question 3:
% Use taylor theorem to get r,v on the first-order:
% r = r +v*delta_t.
% v = v +a*delta_t.
if QRun(3)
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
    PlotFunc(y,z);
    xlabel('Y');
    ylabel('Z');

    %Plot first_order v
    PlotFunc(y_speed,z_speed);
    xlabel('Y Speed');
    ylabel('Z Speed');
end
%% Question 4:
if QRun(4)
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

        taylor_error(t) = sqrt((final_y - t_y(t+1))^2 + (final_z - t_z(t+1))^2);
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
end

%% Question 5,6:

%%%%%%%% Contstants: %%%%%%%%

R = 0.003;                              %Radius in m.
L = 1;                                  %Length in m.
C = 299792458;                          %Speed of light.
m = 1.672621898*(10^(-27));             %Proton mass in MeV.
B = 0.5;                                %Magnetic field.
q = 1.602*(10^(-19));                   %Proton electric change.
w = q*B/m;                              %Max sequence const.
T = abs(2*pi/w);                        %Max time const.
AvgEng = 5;                             %Energy in MeV.
EngDelta = 0.25;                        %Energy delta in MeV.
ProtonMassMeV = 938.272;                %Proton mass in MeV.
E = B*sqrt(2*AvgEng/ProtonMassMeV)*C;   %Electric field according to Q5.

%Acroding to question 5 (in the doc) the fields should be:
if QRun(6)

    PassingEng = AvgEng + 0.1;
    MaxEng = AvgEng + EngDelta;

    %First we need to calculate the speed of 2 energies:
    max_proton_speed = sqrt(2*MaxEng/ProtonMassMeV)*C;
    passing_proton_speed = sqrt(2*PassingEng/ProtonMassMeV)*C;

    delta_t = 0.0000000000001;

    %Set values to zero 
    y = zeros(10000,1);
    z = zeros(10000,1);
    y_speed = zeros(10000,1);
    z_speed = zeros(10000,1);

    % Get initial speed from analitical answer
    y_speed(1) = 0;
    z_speed(1) = max_proton_speed;

    %Repeat the taylor sum for every time jump.
    i = 1;
    while (z(i) < L) && (abs(y(i)) < R) 
        [pos,speed] = TaylorSum(E,B,m,q,delta_t,[0,y(i),z(i)],[0,y_speed(i),z_speed(i)]);
        y(i+1) = pos(2);
        z(i+1) = pos(3);
        y_speed(i+1) = speed(2);
        z_speed(i+1) = speed(3);
        i = i+1;
    end

    %Plot first_order r
    PlotFunc(y,z);
    xlabel('Y');
    ylabel('Z');
end

%% Question 7:
if QRun(7)
    
    MaxEng = AvgEng + EngDelta;
    MinEng = AvgEng - EngDelta;
    avg_speed = sqrt(2*AvgEng/ProtonMassMeV)*C;
    max_speed = sqrt(2*MaxEng/ProtonMassMeV)*C;
    min_speed = sqrt(2*MinEng/ProtonMassMeV)*C;

    passing_dv = [];
    passing_y  = [];

    RepeatNum = 50;
    
    for y0 = -R:(2*R/RepeatNum):R
        for v = min_speed:(max_speed - min_speed)/RepeatNum:max_speed
            
            delta_t = 0.0000000000001;

            %Set values to zero 
            y = zeros(10000,1);
            z = zeros(10000,1);
            y_speed = zeros(10000,1);
            z_speed = zeros(10000,1);

            % Get initial speed from analitical answer
            y(1) = y0;
            z_speed(1) = v;

            i = 1;
            while (z(i) < L) && (abs(y(i)) < R) 
                [pos,speed] = TaylorSum(E,B,m,q,delta_t,[0,y(i),z(i)],[0,y_speed(i),z_speed(i)]);
                y(i+1) = pos(2);
                z(i+1) = pos(3);
                y_speed(i+1) = speed(2);
                z_speed(i+1) = speed(3);
                i = i+1;
            end
            if(z(i) > L || z(i) == L)
                passing_dv = [passing_dv (avg_speed - v)/avg_speed]; 
                passing_y  = [passing_y y0/R];
            end
        end
    end
    %Plot first_order r
    PlotFunc(passing_y,passing_dv);
    xlabel('Y0/R');
    ylabel('deltaV/V0');
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

