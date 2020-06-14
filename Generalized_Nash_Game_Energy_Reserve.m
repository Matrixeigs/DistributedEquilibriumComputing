function result = Generalized_Nash_Game_Energy_Reserve( N )
%GENERALIZED_NASH_GAME Summary of this function goes here
%   Detailed explanation goes here
%Then it the charging plan of each EVs
load Nash_game
SR=[0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 	0.09 ]';
Price_Wholesale_DA=[22.16 	19.41 	19.04 	18.09 	19.06 	22.45 	24.92 	25.41 	26.51 	27.19 	28.28 	30.54 	32.70 	33.04 	35.68 	41.45 	41.97 	44.06 	44.35 	47.66 	40.70 	32.87 	30.16 	26.41 ]';
Price=7*M_t*Price_Wholesale_DA/1000;
SR=7*M_t*SR/1000;
Pr = 3;
%The charging
c_ev = 1.5*max(Price)*ones(N,1);
d_ev = ones(N,1)/1000;
%The decision model of each EV
Lb = zeros(N*T*3,1);
Ub = Pr*ones(N*T*3,1);
%Formulating the A matrix
%sum(Pev)-(Puev)>=minimal departure curve
%sum(Pev)+Pdev<=arrival curve
A_uev = zeros(T,T*3);
%maintance time

for i=1:T
    A_uev(i,1:i) = 1;
    A_uev(i,T+i) = -0.5;
end

A_dev = zeros(T,T*3);
for i=1:T
    A_dev(i,1:i) = 1;
    A_dev(i,2*T+i) = 0.5;
end

A_power_u = zeros(T,T*3);
A_power_d = zeros(T,T*3);
for i=1:T
    A_power_u(i,i) = 1;
    A_power_u(i,T+i) = -1;
    A_power_d(i,i) = 1;
    A_power_d(i,2*T+i) = 1;
end
A_u= zeros(T,T*3);
%maintance time

for i=1:T
    A_u(i,1:i) = 1;
end

A_d = zeros(T,T*3);
for i=1:T
    A_d(i,1:i) = 1;
end

Aev=[-A_uev;A_dev;-A_power_u;A_power_d;-A_u;A_d];
bev= [ ];
for i=1:N
    bev=[bev;-Departure_curve(i,:)';Arrival_curve(i,:)';zeros(T,1);Pr*ones(T,1);-Departure_curve(i,:)';Arrival_curve(i,:)'];
    c(i,:) = [c_ev(i)*ones(T,1);zeros(T,1);zeros(T,1)];
    Q(i,:) = [d_ev(i)*ones(T,1);zeros(T,1);zeros(T,1)]/2;
end

nx=3*T;
ny=6*T;
Price_sg = zeros(nx*N,1);
for i=1:N
    Price_sg((i-1)*nx+1:i*nx) =  [Price; -SR;-SR];
end

Price_variable = sdpvar(nx*N,1);
F = [ ];
for i=1:N
    F = F+[Price_variable((i-1)*nx+1:(i-1)*nx+T)<=Price*2,
        Price_variable((i-1)*nx+1:(i-1)*nx+T)>=0,
        Price_variable((i-1)*nx+T+1:(i-1)*nx+2*T)<=0,
        Price_variable((i-1)*nx+2*T+1:(i-1)*nx+3*T)<=0];
end
M = 1e6;
%The KKT conditions for each EV
X = sdpvar(N*T*3,1);%for each EV, the profiles are given as follows: Pev,Ruev,Rdev
lam_lb = sdpvar(N*T*3,1);
lam_ub = sdpvar(N*T*3,1);
lam = sdpvar(N*ny,1);

B_lb = binvar(N*T*3,1);
B_ub = binvar(N*T*3,1);
B_lam = binvar(N*ny,1);

F = F+[X-Lb>=0,X-Lb<=B_lb*M,Ub-X>=0,Ub-X<=B_ub*M,lam_lb>=0,lam_ub>=0,lam>=0,lam_lb<=(1-B_lb)*M,lam_ub<=(1-B_ub)*M];
for i=1:N
    F = F+[Aev*X(((i-1)*nx+1:i*nx))<=bev(((i-1)*ny+1:i*ny)),bev(((i-1)*ny+1:i*ny))-Aev*X(((i-1)*nx+1:i*nx))<=B_lam(((i-1)*ny+1:i*ny))*M,
        lam(((i-1)*ny+1:i*ny))<=(1-B_lam(((i-1)*ny+1:i*ny)))*M];
    %The KKT system
    %
    F = F+[Price_variable((i-1)*nx+1:i*nx)-c(i,:)'+2*diag(Q(i,:))*X(((i-1)*nx+1:i*nx))-lam_lb(((i-1)*nx+1:i*nx))+lam_ub(((i-1)*nx+1:i*nx))+Aev'* lam(((i-1)*ny+1:i*ny))==0];
end
%For the decision model of power system
H = sparse([ ],[ ],[ ],T,N*nx);
for i=1:T
    H(i,i:nx:nx*N) = 1;
    H(i,2*T+i:nx:nx*N) = 1;
end
Power_limit = 25;
%The power limitation
F = F+[H*X<=Power_limit];
%The price limitation
for i=1:N
    F = F + [c_ev(i)*ones(T,1)-2*d_ev(i)*X((i-1)*nx+1:(i-1)*nx+T)-Price_variable((i-1)*nx+1:(i-1)*nx+T)<=-Price_variable((i-1)*nx+T+1:(i-1)*nx+2*T)];%The price limitation
end
%The objective function is
Q_price = zeros(N*nx,1);
for i=1:N
    Q_price((i-1)*nx+T+1:i*nx) = -1;
end
goal = Price_sg'*X+Price_variable'*Q_price;
%Deploy the objective functions
option = sdpsettings('solver','cplex');
sol = solvesdp(F,goal,option);

X = double(X);
Price_variable = double(Price_variable);
sg_utility = (Price_variable-Price_sg)'*X;

%The utility of each EV
EV_utility=zeros(N,1);
for i=1:N
    EV_utility(i)=c_ev(i)*sum(X((i-1)*nx+1:(i-1)*nx+T))-X((i-1)*nx+1:(i-1)*nx+T)'*diag(Q(i,1:T))*X((i-1)*nx+1:(i-1)*nx+T)-Price_variable((i-1)*nx+1:i*nx)'*X((i-1)*nx+1:i*nx);
end
result=zeros(N,nx);
for i=1:N
    result(i,:) = X((i-1)*nx+1:i*nx);
end
Energy_plan=result(:,1:T);
Ru=result(:,T+1:2*T);
Rd=result(:,2*T+1:3*T);
%By using the gurobi model
Energy_price = zeros(N,T);
Ru_price = zeros(N,T);
Rd_price = zeros(N,T);
for i=1:N
    Energy_price(i,:)=Price_variable((i-1)*nx+1:(i-1)*nx+T);
    Ru_price(i,:)=-Price_variable((i-1)*nx+T+1:(i-1)*nx+2*T);
    Rd_price(i,:)=-Price_variable((i-1)*nx+2*T+1:(i-1)*nx+3*T);
end

social_welfare= sg_utility+sum(EV_utility);

save Result

end

