function result =Nash_Game( N )
%GENERALIZED_NASH_GAME Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    N=100;%There are 1000 EVs
end

%Generalized Nash Game Approach
%The competition among different EV users are converted into a potential
%game, i.e., a unique NE exits.
%1) Generate the integration time of each EVs
%The arrivial time
%The energy required
%The departure time
x1=0:0.01:7.5;
y1=1/1.2/pi./(1+5*(x1-7.3).^2);
x2=7.51:0.01:23.99;
y2=1/1.05/pi./(1+6*(x2-8.2).^2);
x=[x1,x2];
y=[y1,y2];
y_m_s=100*y/sum(y);%???
y_e_e=normpdf(x,18.2,0.8);
for i=1:24
    Y_s(i)=sum(y_m_s((i-1)*100+1:i*100));
    Y_e(i)=sum(y_e_e((i-1)*100+1:i*100));
end
clear x y;%????
Y_s=Y_s/100;%????
Y_e=Y_e/100;%????
%Generate the energy demand of each EV

T_shift=15;
T=24;
M_t=zeros(T,T);
for i=1:T
    if i<=T-T_shift
        M_t(i,i+T_shift)=1;
    else
        M_t(i,i-T+T_shift)=1;
    end
end

Y_s=M_t*Y_s';
Y_e=M_t*Y_e';

Y_s_cum=cumsum(Y_s);
Y_e_cum=cumsum(Y_e);

Cr=30;%????
Pr=3;%??????
eff_c=0.9;%????
Delta_t=1;
%Generate the arrival curve and departure curve of each EV
%Generate each EV's information:arrival time, departure time, and energy
%required
index=1:T;
for i=1:N
    Arrival_slot = max(index(rand>Y_e_cum));
    if ~isempty(Arrival_slot)
        EV(i).Arrival_time=Arrival_slot;
    else
        EV(i).Arrival_time=1;
    end
    Departure_slot = max(index(rand>Y_s_cum));
    if ~isempty(Departure_slot)
        EV(i).Departure_time=Departure_slot;
    else
        EV(i).Departure_time=1;
    end
    EV(i).Energy_required=2*gamrnd(5,3)*0.15/eff_c;
    if EV(i).Departure_time<= EV(i).Arrival_time
        EV(i).Departure_time= EV(i).Arrival_time;
    end
    if EV(i).Energy_required>=(EV(i).Departure_time-EV(i).Arrival_time)*Pr*eff_c*Delta_t
        EV(i).Energy_required = (EV(i).Departure_time-EV(i).Arrival_time)*Pr*eff_c*Delta_t;
        if EV(i).Energy_required == 0
            EV(i).Energy_required = rand;%A small random value is assigned to this field.
        end
    end
    
    EV(i).minimal_charging_duration = ceil((EV(i).Energy_required)/eff_c/Pr/Delta_t);
    
    Profile = zeros(1,EV(i).minimal_charging_duration);
    E0 = EV(i).Energy_required;
    j = 1;
    while E0>0 && j<=EV(i).minimal_charging_duration
        Profile(j) = min(E0/eff_c/Delta_t, Pr);
        E0 = E0 -Profile(j) *eff_c*Delta_t;
        j = j + 1;
    end
    
    EV(i).maximal_charging_profile = Profile;
    
end
Arrival_curve=zeros(N,T);
Departure_curve=zeros(N,T);
for i=1:N
    Arrival_curve(i, 1:EV(i).Arrival_time-1) = 0;
    Arrival_curve(i, EV(i).Arrival_time:EV(i).Arrival_time+EV(i).minimal_charging_duration-1) = cumsum(EV(i).maximal_charging_profile);
    Arrival_curve(i, EV(i).Arrival_time+EV(i).minimal_charging_duration : T) = sum(EV(i).maximal_charging_profile);
    Departure_curve(i, 1:EV(i).Departure_time-EV(i).minimal_charging_duration) = 0;
    Departure_curve(i, EV(i).Departure_time-EV(i).minimal_charging_duration+1:EV(i).Departure_time) = cumsum(sort(EV(i).maximal_charging_profile,'ascend'));
    Departure_curve(i, EV(i).Departure_time+1:T) =  sum(EV(i).maximal_charging_profile);
end
%Then it the charging plan of each EVs
Price=[0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;1.22830000000000;1.22830000000000;1.22830000000000;1.22830000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;1.33770000000000;1.33770000000000;1.33770000000000;1.33770000000000;0.778500000000000;0.778500000000000;0.353900000000000];
Price=M_t*Price;
save EV
%Deploy the objective functions
Price_b=ones(T,1)*0.01;
%Price_b=zeros(T,1)*0.01;
H = sparse([ ],[ ],[ ],T,N*T);
for i=1:T
    H(i,i:T:N*T) = 1;
end

%By using the gurobi model
Lb=zeros(N*T,1);
Ub=zeros(N*T,1);
for i=1:N
    Ub((i-1)*T+EV(i).Arrival_time:(i-1)*T+EV(i).Departure_time)=Pr;
end
A = [ ];
bu = [ ];
bl = [ ];
B = [ ];
N_lam = zeros(N,1);
for i=1:N
    A(i).Cum = [ ];
    for j=EV(i).Arrival_time:EV(i).Departure_time
        A(i).Cum = [A(i).Cum ; sparse(ones(1,j-EV(i).Arrival_time+1),EV(i).Arrival_time:j,ones(1,j-EV(i).Arrival_time+1),1,T)];
    end
    A(i).Cum = [ A(i).Cum;eye(T,T);-A(i).Cum;-eye(T,T);];
    B(i).b=[Arrival_curve(i,EV(i).Arrival_time:EV(i).Departure_time)';Ub((i-1)*T+1:i*T);-Departure_curve(i,EV(i).Arrival_time:EV(i).Departure_time)';-Lb((i-1)*T+1:i*T);];
    N_lam(i)=length(B(i).b);
end

% M=1e3;
N_lam_sum=sum(N_lam);
% X = sdpvar(N*T,1);
N_lam=cumsum(N_lam);
% % 
% lam_b=sdpvar(N_lam_sum,1);
% 
% B_lam_b=binvar(N_lam_sum,1);
% F= [ ];
% for i = 1 : N
%     %The objective function for each player
%     F = F + [ A(i).Cum * X((i-1)*T+1:i*T) - B(i).b<=0 ];
%     
%     if i == 1
%         F = F + [ B(i).b - A(i).Cum * X((i-1)*T+1:i*T) <= B_lam_b(1:N_lam(1)) * M];
%         F = F + [ lam_b(1:N_lam(1)) <=(1- B_lam_b(1:N_lam(1))) * M];
%         F = F + [ lam_b(1:N_lam(1)) >= 0];
%     else
%         F = F + [ B(i).b - A(i).Cum * X((i-1)*T+1:i*T)  <= B_lam_b(N_lam(i-1)+1:N_lam(i))*M];
%         F = F + [ lam_b(N_lam(i-1)+1:N_lam(i)) <=(1- B_lam_b(N_lam(i-1)+1:N_lam(i))) * M];
%         F = F + [ lam_b(N_lam(i-1)+1:N_lam(i)) >= 0];
%     end
%     %The KKT conditions
%     if i==1
%         F= F + [Price+diag(Price_b/2)*H*X + diag(Price_b/2)*X((i-1)*T+1:i*T) + A(i).Cum'*lam_b(1:N_lam(1)) == 0];
%     else
%         F= F + [Price+diag(Price_b/2)*H*X + diag(Price_b/2)*X((i-1)*T+1:i*T) + A(i).Cum'*lam_b(N_lam(i-1)+1:N_lam(i)) == 0];
%     end
%     
% end
% option = sdpsettings('solver','gurobi');
% sol = solvesdp(F,[ ],option);
% % 
% X = double(X);
% X_sum =H*X;
% Price_ng=Price+diag(Price_b)*X_sum;
% plot(Price_ng);
%Results analysis
% plot(Price_ng)
% hold on
% plot(diag(Price_b/2)*X_sum,'r');
% figure
% plot(X_sum/1000);
% hold on
% plot(Price)
% hold on
% plot(Price0,'r')
A_0 = sparse([ ],[ ],[ ], N_lam_sum,N*T);
bu = [ ];
for i=1:N
    if i==1
         A_0(1:N_lam(i),(i-1)*T+1:i*T)= A(i).Cum;
    else
        A_0(N_lam(i-1)+1:N_lam(i),(i-1)*T+1:i*T)= A(i).Cum;
    end
    bu = [bu;B(i).b];
end

%The distributed algorithms
gap =inf;
Time_cost=0;
params.OutputFlag=0;
clear models
model.obj = H'*Price;
model.Q = sparse(H'*diag(Price_b/4)*H)+sparse(1:N*T,1:N*T,repmat(Price_b/4,N,1),N*T,N*T);
model.modelsense = 'min';
model.A = A_0;
model.rhs = bu;
model.sense = repmat('<',length(bu),1);
%The constrain set

results  = gurobi(model,params);
model_social=model;
model_social.Q = sparse(H'*diag(Price_b/2)*H);
results_social  = gurobi(model_social,params);


X0=results.x;

X0_sum=H*X0;

Price_ng = Price + diag(Price_b)*X0_sum;

social_welfare_0=Price'*X0_sum+X0_sum'*diag(Price_b/2)*X0_sum;
% social_welfare=Price'*X_sum+X_sum'*diag(Price_b/2)*X_sum;
Departure_curve_real = zeros(N,T);

%The departure curve depiction
for i=1:N
    Departure_curve_real(i,:)=cumsum(X0((i-1)*T+1:i*T));
end
Uncontrolled_charging = zeros(N,T);
for i=1:N
    Uncontrolled_charging(i, EV(i).Arrival_time:EV(i).Arrival_time+EV(i).minimal_charging_duration-1)=EV(i).maximal_charging_profile;
end

%The social welfare
X_social=results_social.x;

X_social_sum=H*X_social;
social_welfare=Price'*X_social_sum+X_social_sum'*diag(Price_b/2)*X_social_sum;
Price_welfare=Price + diag(Price_b)*X_social_sum;

% gap = norm(X0_sum-X_sum);
save Nash_game

end