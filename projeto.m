% 1.5
% Define symbolic variables
pkg load symbolic


syms phi theta psi real
syms dphi dtheta dpsi real
syms omega_x omega_y omega_z real
syms velocity_x velocity_y velocity_z real

%Define constants
beta = 0.1021;
m = 1539;
g = 3.71;

% Set angle values
phi = 0;
theta = 0;
psi = 0;



rotx = [	1	,	0			,	 0
			0	,	cos(phi)	,	-sin(phi)
			0	,	sin(phi)	,	 cos(phi)	];

roty = [	 cos(theta)	,	0	,	sin(theta)
			0				,	1	,	0
			-sin(theta)	,	0	,	cos(theta)	];

rotz = [	cos(psi)	,	-sin(psi)	,	0
			sin(psi)	,	 cos(psi)	,	0
			0			,	 0				,	1	];

% First row elements  1.2.1

% Define velocities(v)
velocity_x = 0 ;
velocity_y = 0 ;
velocity_z = 89;

v=[velocity_x;velocity_y;velocity_z]

% Define the matrix df_p/dv
dfp_dv = rotz * roty * rotx;

% Define the matrix df_p/dphi
dfp_dphi= [0,     0    ,    0     ;
           0, -sin(phi), -cos(phi);
           0, cos(phi) , -sin(phi)] * roty * rotx;

% Define the matrix df_p/dtheta
dfp_dtheta = rotx * [-sin(theta),     0    ,    cos(theta);
                          0     ,     0    ,      0       ;
                     -cos(theta),     0    , -sin(theta)] * rotz;

% Define the matrix df_p/dpsi
dfp_dpsi = rotx * roty * [-sin(psi), -cos(psi)    ,       0     ;
                           cos(psi), -sin(psi)    ,       0     ;
                              0    ,     0        ,       0     ];

% Compute the overall derivative df_p/dlambda
dfp_dlambda = [dfp_dphi; dfp_dtheta; dfp_dpsi] * [velocity_x, velocity_y, velocity_z]';

% Display the results
disp('df_p/dphi:');
disp(dfp_dphi);

disp('df_p/dtheta:');
disp(dfp_dtheta);

disp('df_p/dpsi:');
disp(dfp_dpsi);

disp('df_p/dv:');
disp(dfp_dv);

disp('df_p/dlambda:');
disp(dfp_dlambda);



% Second row elements 1.2.2

% Define angular velocities(w)
omega_x= 0 ;
omega_y= 0 ;
omega_z= 0 ;
omega=[omega_x;omega_y;omega_z]
% Define the matrix df_v/dv
dfv_dv = -skew(omega) - (2 * beta) / m;

% Define the matrix df_v/dlambda
%z_I
z_I = [0; 0; 1];

dfv_dlambda = g * rotx * roty * [-sin(psi), -cos(psi), 0;
                                  cos(psi), -sin(psi), 0;
                                      0   ,     0    , 0];

% Define the matrix df_v/domega
dfv_domega = skew(v);

% Display the results
disp('df_v/dv:');
disp(dfv_dv);

disp('df_v/dlambda:');
disp(dfv_dlambda);

disp('df_v/domega:');
disp(dfv_domega);



% Third row elements  1.2.3



% Define the matrix df_lambda/dphi
dflambda_dphi = [0,       cos(phi)*tan(theta)        ,      -sin(phi)*tan(theta)         ;
                 0,           -sin(phi)              ,            -cos(phi)              ;
                 0, cos(theta)*cos(phi)/cos(theta).^2, -cos(theta)*sin(phi)/cos(theta).^2];

% Define the matrix df_lambda/dtheta
dflambda_dtheta = [0,       cos(phi)*sec(theta).^2     , cos(phi)*sec(theta).^2        ;
                   0,                  0               ,               0               ;
                   0, sin(phi)*sec(theta).^2*sin(theta), cos(phi)*sec(theta)*tan(theta)];

% Define the matrix df_lambda/dpsi
dflambda_dpsi = zeros(3, 3);

% Define the matrix df_lambda/domega
dflambda_domega = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
                   0,      cos(phi)      ,     -sin(phi)      ;
                   0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

% Compute the overall derivative df_lambda/dlambda
dflambda_dlambda = [dflambda_dphi, dflambda_dtheta, dflambda_dpsi]' * [omega_x; omega_y; omega_z];

% Display the results
disp('df_lambda/dphi:');
disp(dflambda_dphi);

disp('df_lambda/dtheta:');
disp(dflambda_dtheta);

disp('df_lambda/dpsi:');
disp(dflambda_dpsi);

disp('df_lambda/domega:');
disp(dflambda_domega);

disp('df_lambda/dlambda:');
disp(dflambda_dlambda);



% Fourth row elements  1.2.3

% Define J
Jxx = 2332;
Jyy = 2551;
Jzz = 2090;

% Define the matrix df_omega/domega
dfomega_domega = [0, (omega_z * (Jxx - Jyy)) / Jxx, (omega_y * (Jzz - Jyy)) / Jxx;
                  (omega_z * (Jxx - Jzz)) / Jyy, 0, (omega_x * (Jxx - Jzz)) / Jyy;
                  (omega_y * (Jyy - Jxx)) / Jzz, (omega_z * (Jxx - Jzz)) / Jyy, 0];

% Display the results
disp('df_omega/domega:')
disp(dfomega_domega);




