%state machine

clear
m = 1539;
m2b = m/3;
zI=[0;0;1]

syms phi theta psi real


A1a = [
  zeros(3), eye(3), [0, -89, 0; 89, 0, 0; 0, 0, 0], zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; 3.71, 0, 0; 0, 0, 0], [0, -89, 0; 89, 0, 0; 0, 0, 0];
  zeros(3), zeros(3), zeros(3), eye(3);
  zeros(3), zeros(3), zeros(3), zeros(3)
];

A1b = [
  zeros(3), eye(3), [0, -45, 0; 45, 0, 0; 0, 0, 0], zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; 3.71, 0, 0; 0, 0, 0], [0, -89, 0; 89, 0, 0; 0, 0, 0];
  zeros(3), zeros(3), zeros(3), eye(3);
  zeros(3), zeros(3), zeros(3), zeros(3)
];

A2a = [
  zeros(3), eye(3), zeros(3), zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; 3.71, 0, 0; 0, 0, 0], zeros(3);
  zeros(3), zeros(3), zeros(3), eye(3);
  zeros(3), zeros(3), zeros(3), zeros(3)
];

syms phi theta psi real
% Set angle values
phi = 90;
theta = 0;
psi = 0;


A2b = [
  zeros(3), [1, 0, 0; 0, 0, -1; 0, 1, 0], [0, 3.2, -2.4; 0, 0, 0; -4, 0, 0], zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; -1.6624, 0, 0; 3.3167, 0, 0], [0, 0, 4; 0, 0, 0; -4, 0, 0];
  zeros(3), zeros(3), zeros(3), [1, 0, 0; 0, -0.448, -0.894; 0, 0.894, 0.4481];
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

B2b = [
    zeros(3, 4);
    (1/m2b) * [
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

% get LQR-1a controller
Q1a = blkdiag(9*eye(3),10*eye(3),20*eye(3),30*eye(3));
R1a = 0.01*eye(4);
Klqr1a = lqr(A1a,B,Q1a,R1a);
lbd_CL_lqr1a = eig(A1a-B*Klqr1a); 
if any(real(lbd_CL_lqr1a) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

B1 = zeros(12,12);
B2 = B;
W1_1a = sqrt(Q1a);
W2_1a = sqrt(R1a);
C1_1a = [W1_1a];
D12_1a = [zeros(8,4); W2_1a];
D11 = zeros(12,12);
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);
B0 = [B1, B2];
C01a = [C1_1a; C2];
D01a = [D11, D12_1a; D21, D22];


% get LQR-1b controller
Q1b = blkdiag(9*eye(3),10*eye(3),20*eye(3),30*eye(3));
R1b = 0.01*eye(4);
Klqr1b = lqr(A1b,B,Q1b,R1b);
lbd_CL_lqr1b = eig(A1b-B*Klqr1b); 
if any(real(lbd_CL_lqr1b) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

B1 = zeros(12,12);
B2 = B;
W1_1b = sqrt(Q1b);
W2_1b = sqrt(R1b);
C1_1b = [W1_1b];
D12_1b = [zeros(8,4); W2_1b];
D11 = zeros(12,12);
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);
B0 = [B1, B2];
C01b = [C1_1b; C2];
D01b = [D11, D12_1b; D21, D22];

% get LQR-2a controller
Q2a = blkdiag(9*eye(3),10*eye(3),20*eye(3),30*eye(3));
R2a = 0.01*eye(4);
Klqr2a = lqr(A2a,B,Q2a,R2a);
lbd_CL_lqr2a = eig(A2a-B*Klqr2a); 
if any(real(lbd_CL_lqr2a) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

B1 = zeros(12,12);
B2 = B;
W1_2a = sqrt(Q2a);
W2_2a = sqrt(R2a);
C1_2a = [W1_2a];
D12_2a = [zeros(8,4); W2_2a];
D11 = zeros(12,12);
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);
B0 = [B1, B2];
C02a = [C1_2a; C2];
D02a = [D11, D12_2a; D21, D22];

% get LQR-2b controller
Q2b = blkdiag(9*eye(3),10*eye(3),20*eye(3),30*eye(3));
R2b = 0.01*eye(4);
Klqr_2b = lqr(A2b,B2b,Q2b,R2b);
lbd_CL_lqr2b = eig(A2b-B2b*Klqr_2b); 
if any(real(lbd_CL_lqr2b) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

B1 = zeros(12,12);
W1_2b = sqrt(Q2b);
W2_2b = sqrt(R2b);
C1_2b = [W1_2b];
D12_2b = [zeros(8,4); W2_2b];
D11 = zeros(12,12);
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);
B02b = [B1, B2b];
C02b = [C1_2b; C2];
D02b = [D11, D12_2b; D21, D22];


%P matrix %
P1a = ss(A1a, B0, C01a, D01a);
P1b = ss(A1b, B0, C01b, D01b);
P2a = ss(A2a, B0, C02a, D02a);
P2b = ss(A2b, B02b, C02b, D02b);


nmeas = 12; % Number of measured outputs
ncont = 4; % Number of control inputs

[Kinf1a,CLinf1a,gammainf1a,info_inf1a] = hinfsyn(P1a,nmeas,ncont);
poles_CLinf1a = pole(CLinf1a);
if any(real(poles_CLinf1a) >= 0), disp('CL-1a system with Hinf controller not stable'); else, disp('CL-1a system with Hinf controller is stable'); end

[Kinf1b,CLinf1b,gammainf1b,info_inf1b] = hinfsyn(P1b,nmeas,ncont);
poles_CLinf1b = pole(CLinf1b);
if any(real(poles_CLinf1b) >= 0), disp('CL-1b system with Hinf controller not stable'); else, disp('CL-1b system with Hinf controller is stable'); end

[Kinf2a,CLinf2a,gammainf2a,info_inf2a] = hinfsyn(P2a,nmeas,ncont);
poles_CLinf2a = pole(CLinf2a);
if any(real(poles_CLinf2a) >= 0), disp('CL-2a system with Hinf controller not stable'); else, disp('CL-2a system with Hinf controller is stable'); end

[Kinf2b,CLinf2b,gammainf2b,info_inf2b] = hinfsyn(P2b,nmeas,ncont);
poles_CLinf2b = pole(CLinf2b);
if any(real(poles_CLinf2b) >= 0), disp('CL-2b system with Hinf controller not stable'); else, disp('CL-2b system with Hinf controller is stable'); end



% Hinf controller
Dt = 0.05;
t = 0:Dt:100;
r = [0;0;-500;0;0;45;0;0;0;0;0;0]*(t>=0);
r1a = [0;0;-500;0;0;45;0;0;0;0;0;0]*(t>=10);
r1b = [0;0;-22;0;0;0;0;0;0;0;0;0]*(t>=10);
r2a = [0;0;-22;0;4;0;90;0;0;0;0;0]*(t>=10);
r2b = [0;100;-22;0;4;0;0;90;0;0;0;0]*(t>=10); %onde queremos que ele vá (referencia) OP1b
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-2100;0;0;89;0;0;0;0;0;0]
C=eye(12)
for k = 1:NSim
    
    y(:,k) = C*x(:,k);
    z=y(3,:);
    phi=y(7,:);
    aux=0;
    
    if z<-500
        x(:,1) = [0;0;-2100;0;0;89;0;0;0;0;0;0];
        r(:,k) = [0;0;-500;0;0;45;0;0;0;0;0;0]*(t(k)>=0);
        Kinf=Kinf1a;
        A=A1a;
        aux=1;
    elseif  ((z>=-500) & (z~=-22) & (aux==1))
        x(:,1) = [0;0;-500;0;0;45;0;0;0;0;0;0];
        r(:,k) = [0;0;-22;0;0;0;0;0;0;0;0;0]*(t(k)>=0);
        Kinf=Kinf1b;
        A=A1b;
        aux=2;
    elseif ((z==-22) & (phi~=0) & (aux==2))
         Kinf=Kinf2b;
         x(:,1) = [0;0;-22;0;4;0;0;90;0;0;0;0]; %onde ele está condicoes iniciais OP1a
         r(:,k) = [0;100;-22;0;4;0;0;90;0;0;0;0]*(t(k)>=0);
         A=A2b;
         B=B2b;
    else
        x(:,1) = [0;0;-22;0;0;0;0;0;0;0;0;0];
        r(:,k) = [0;0;-22;0;4;0;90;0;0;0;0;0]*(t(k)>=0);
        Kinf=Kinf2a;
        A=A2a;
    end
    
    
    % get measurements:

    %y1(:,k) = C*x(:,k);
    
    % get control action:

    v = -[(r(:,k)-y(:,k))];
    u(:,k) = Kinf.C*v; % approximation for LQR-like performance of Hinf

    % simulate system:
    x_dot = A*x(:,k) + B*u(:,k); % system derivatives
    xp = x(:,k) + Dt*x_dot; % integrate system state
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(9041002);
plot(t,u);
grid on;
xlabel('t[s]');
ylabel('u(t)');

figure(9041003);
plot(t,y(3,:),t,r(3,:)); %referencia em z = 3
grid on;
xlabel('t[s]');
legend('z','r_z');

figure(9041004);
plot(t,y(7,:),t,r(7,:)); %referencia em z = 3
grid on;
xlabel('t[s]');
legend('\phi','r_{\phi}');

figure(9041005);
plot(t,y(2,:)); %referencia em z = 3
grid on;
xlabel('t[s]');
legend('y');






