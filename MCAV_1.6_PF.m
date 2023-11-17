%% Controlability and Observability
pkg load control

m = 1539;
zI=[0;0;1]

A = [
  zeros(3), eye(3), [0, -89, 0; 89, 0, 0; 0, 0, 0], zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; 3.71, 0, 0; 0, 0, 0], [0, -89, 0; 89, 0, 0; 0, 0, 0];
  zeros(3), zeros(3), zeros(3), eye(3);
  zeros(3), zeros(3), zeros(3), zeros(3)
];

B = [
    zeros(3, 4);
    (1/m) * [
        [0; 0.342; -0.940], [0; 0.342; -0.940], [0; -0.342; -0.940], [0; -0.342; -0.940]
    ];
    zeros(3, 4);
    [
        [6.63e-4; -4.97e-4; 2.21e-4], [6.63e-4; 4.97e-4; 2.21e-4], [-6.63e-4; 4.97e-4; -2.21e-4], [-6.63e-4; -4.97e-4; 2.21e-4]
    ]
];
C = [   eye(3)    , zeros(3)   , zeros(3) , zeros(3)
        zeros(1,3), zeros(1,3) , zI'      , zeros(1,3)  ];
D = zeros(4);

sys = ss(A,B,C,D);

[Vj,Jor] = schur(A),
% [Vj,Jor] = jordan(sym(A)),
[V,DL,W] = eig(A),
mode_obs = C*V,
mode_ctrl = W'*B,


