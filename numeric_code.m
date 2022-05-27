%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 27/05/2022   %
%  Electricity Numeric Ex   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Loaded particles influenced by electric and magnetic fields

%%%%%%%% Contstants: %%%%%%%%

E = 1;              %Electric field const.
B = -1;             %Magnetic field const.
q = 1;              %particle charge const.
m = 1;              %particle mass const.
w = q*B/m;          %Max sequence const.
T = 2*pi/w;         %Max time const.

%%%%%%%% Questions: %%%%%%%%
%% Question 2:

%Set initial values and consts:
t = linspace(0,T);

% The particle's movement vs time goes as followes (explained in doc):
% x = 0;
% y = (2E/wb)*cos(wt) - (2E/wb).
% z = (2E/wb)*sin(wt) + (Et/b).

y = (2*E/(w*B))*cos(w.*t) - 2*E/(w*B);
z = (2*E/(w*B))*sin(w.*t) - E.*t./B;

PlotFunc(y,z);
xlabel('Y');
ylabel('Z');

%% Question 3:


%%%%%%%% Functions: %%%%%%%%

function p = PlotFunc(x,y)
    figure
    hold on
    box on
    grid
    plot(x,y)
    hold off
end 

