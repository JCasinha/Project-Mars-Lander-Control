%non linear model (1.4)

pkg load symbolic

%syms phi theta psi real
%syms omega_x omega_y omega_z real
%syms velocity_x velocity_y velocity_z real

%define constants
g = 3.71;  % Gravitational constant on Mars (m/s^2)
Cp=2.1
rho=0.0120
Az=8.1
Ay=5.94
bz = 0.5*Cp*rho*Az;
by = 0.5*Cp*rho*Ay;      % Drag coefficient
m = 513;      % Mass of the lander (kg) 1/2 of the rover


%nonlinear parameters
phi=0;
theta=0;
psi=0;

vx=0;
vy=0;
vz=89; %m/s

w1=0;
w2=0;
w3=0;


%rotation matrix

rotx = [	1	,	0			,	 0
			0	,	cos(phi)	,	-sin(phi)
			0	,	sin(phi)	,	 cos(phi)	];

roty = [	 cos(theta)	,	0	,	sin(theta)
			0				,	1	,	0
			-sin(theta)	,	0	,	cos(theta)	];

rotz = [	cos(psi)	,	-sin(psi)	,	0
			sin(psi)	,	 cos(psi)	,	0
			0			,	 0				,	1	];

R = rotx*roty*rotz;

%Q matrix

Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
     0, cos(phi), -sin(phi);
     0, sin(phi)/cos(theta), cos(phi)/cos(theta);]


%J matrix

J=[2332,0,0;0,2551,0;0,0,2090]

%non linear vectors
lbd=[phi;theta;psi]
v=[vx;vy;vz]
om=[w1;w2;w3]
zI=[0;0;1]

%C matrix

C = [   eye(3)    , zeros(3)   , zeros(3) , zeros(3)
        zeros(1,3), zeros(1,3) , zI'      , zeros(1,3)  ];

% simulation parameters
Dt = 0.01;
t = 0:Dt:10;
nx = 12;
ny = 4;
x0 = zeros(nx,1);


% simulate nonlinear system
u_L = [0.2*m*g;0.01;0.01;0.01]*(t>=0);
Te = m*g; % equilibirum thrust
u_NL = [Te;0;0;0]*ones(size(t)) + u_L;
Nsim = length(t);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);
x(:,1) = x0;


for k = 1:Nsim
    % prepare variables:
    T = u_NL(1,k);
    np = u_NL(2:4,k);

    % compute state derivative:
    p_dot = R*v;
    v_dot = -skew(om)*v + g*R'*zI - 1/m*T;
    lbd_dot = Q*om;
    om_dot = -inv(J)*skew(om)*J*om + inv(J)*np;
    x_dot = [p_dot;v_dot;lbd_dot;om_dot];

    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot;

    % compute current output:
    y(:,k) = C*x(:,k);
    asd = 1;
end

figure(90322);
plot(t, y(3:4, :), '-.'); % Nonlinear system only
grid on;
xlabel('Time (s)');
ylabel('Z Position (m)');
legend('Lander Velocity', 'Resting State');
title('Nonlinear Model Simulation');



