function result=Bilevel_final_fixed(Emax_interval,Emin_interval,T,Peak_load,NPro,N,Pr,W)
tic;
%The methods are not fully demonstrated!
Price=[0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;0.353900000000000;1.22830000000000;1.22830000000000;1.22830000000000;1.22830000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;0.778500000000000;1.33770000000000;1.33770000000000;1.33770000000000;1.33770000000000;0.778500000000000;0.778500000000000;0.353900000000000];
T_shift=15;%默认起始控制时段为第30个时间点
if isempty(T_shift)
    T_shift=15;
end
M_t=zeros(T,T);
for i=1:T
    if i<=T-T_shift
        M_t(i,i+T_shift)=1;
    else
        M_t(i,i-T+T_shift)=1;
    end
end
Price=M_t*Price;

Emax=cumsum(Emax_interval*M_t');%应该先平移后相加
Emin=cumsum(Emin_interval*M_t');%各时段可调电量上下界确定

%%计算模块.
%分配到各个工区
%EV比例
NA=5;%下属子区域数量
%Price=repmat(Price,NA,1);

load Data;
Loadprofile=Pd*Peak_load;%各时段的负荷分布
Loadprofile=Loadprofile*M_t';
Loadprofile=Loadprofile';
%充电设施数量
Charger_pro=[0.1,0.1,0.2,0.3,0.3]*NPro;
%电动汽车比例
EV_pro=[0.1,0.1,0.2,0.3,0.3];%各个区域内电动汽车的数量
Pmax_area=reshape(ones(T,1)*Pr*Charger_pro*N,NA*T,1);

Emax_area=reshape(Emax'*EV_pro,T*NA,1);
Emin_area=reshape(Emin'*EV_pro,T*NA,1);
Emax=Emax';
Emin=Emin';
EV_un=zeros(T,NA);
for i=1:NA
    temp=zeros(T,1);
    Emax_temp=[Emax_area((i-1)*T+1);diff(Emax_area((i-1)*T+1:i*T))];%最大累计电量
    Pmax_temp=Pmax_area((i-1)*T+1:i*T);
    for j=1:T
        if j==1
            temp(j)=min(Emax_temp(j),Pmax_temp(j));
        else
            Emax_temp(j)=Emax_temp(j)+Emax_temp(j-1)-temp(j-1);
            temp(j)=min(Emax_temp(j),Pmax_temp(j));
        end
    end
    EV_un(:,i)=temp;
end

%约束条件生成
Loadpro=[0.4 0.2 0.15 0.15 0.1];
Load_area=reshape(Loadprofile*Loadpro,T,NA);%各区域的空间负荷分布
load PA
Load_area=PA*Peak_load;
Loadprofile=sum(Load_area);
Loadprofile=Loadprofile';
Loadprofile=M_t*Loadprofile;
Load_area=M_t*Load_area';
%Upper Level Control
x=sdpvar(T,1);%各时段的充电功率，已经是充电功率，不需要折合
y=sdpvar(T*NA,1);%为各区域的充电功率，也已经折合
%充电功率的累积关系
F0=[x>=0,x<=N*Pr];
A=tril(ones(T,T));%累积充电功率矩阵
F0=[F0,A*x<=Emax,A*x>=Emin];
H=eye(T,T)-1/T*ones(T,T);%二次矩阵；
goal0=x'*H*x+2*Loadprofile'*H*x;
%goal0=Price'*x;
%re=solvesdp(F0,goal0);

Atrans=zeros(T,T*NA);
for i=1:T
    Atrans(i,i:T:end)=1;
end
Acum=zeros(NA*T,NA*T);
for i=1:NA
    Acum((i-1)*T+1:i*T,(i-1)*T+1:i*T)=A;
end

F1=[y>=0,y<=Pmax_area,Acum*y>=Emin_area,Acum*y<=Emax_area];
y1=sdpvar(NA*T,1);
%权重的设定影响较大
u1=sdpvar(size(y,1),1);
M1=binvar(size(y,1),1);
u2=sdpvar(size(y,1),1);
M2=binvar(size(y,1),1);
u3=sdpvar(size(Acum,1),1);
M3=binvar(size(Acum,1),1);
u4=sdpvar(size(Acum,1),1);
M4=binvar(size(Acum,1),1);
M0=100000;
F1=[F1,u1>=0,u2>=0,u3>=0,u4>=0];%拉格朗日因子
%The problem within each area
W1=0;
for i=1:NA
    %Lagrange Equation
    if i>NA
         F1=[F1,2*A'*A*y((i-1)*T+1:i*T)-2*A'*Emax_area((i-1)*T+1:i*T)+W1*2*y((i-1)*T+1:i*T)-2*W1*y1((i-1)*T+1:i*T)-u1((i-1)*T+1:i*T)+u2((i-1)*T+1:i*T)-A'*u3((i-1)*T+1:i*T)+A'*u4((i-1)*T+1:i*T)==0];%Dos
    else
        F1=[F1,Price+W1*2*y((i-1)*T+1:i*T)-2*W1*y1((i-1)*T+1:i*T)-u1((i-1)*T+1:i*T)-u1((i-1)*T+1:i*T)+u2((i-1)*T+1:i*T)-A'*u3((i-1)*T+1:i*T)+A'*u4((i-1)*T+1:i*T)==0];%Cost
    end
    F1=[F1,u1((i-1)*T+1:i*T)<=M0*M1((i-1)*T+1:i*T),y((i-1)*T+1:i*T)<=M0*(1-M1((i-1)*T+1:i*T))];
    F1=[F1,u2((i-1)*T+1:i*T)<=M0*M2((i-1)*T+1:i*T),Pmax_area((i-1)*T+1:i*T)-y((i-1)*T+1:i*T)<=M0*(1-M2((i-1)*T+1:i*T))];
    F1=[F1,u3((i-1)*T+1:i*T)<=M0*M3((i-1)*T+1:i*T),A*y((i-1)*T+1:i*T)-Emin_area((i-1)*T+1:i*T)<=M0*(1-M3((i-1)*T+1:i*T))];
    F1=[F1,u4((i-1)*T+1:i*T)<=M0*M4((i-1)*T+1:i*T),Emax_area((i-1)*T+1:i*T)-A*y((i-1)*T+1:i*T)<=M0*(1-M4((i-1)*T+1:i*T))];
end


%x=Atrans*y;
y1=Atrans*y;
%W=0.01;
goal1=goal0+W*(y1'*y1-2*x'*y1+x'*x);

option=sdpsettings('solver','cplex','verbose',0);
re=solvesdp(F0+F1,goal1,option);

u1=double(u1);
M2=double(M1);
u2=double(u2);
M2=double(M2);
u3=double(u3);
M3=double(M3);
u4=double(u4);
M4=double(M4);
x=double(x);
y=double(y);
y1=double(y1);
result.Cost=Price'*y1;
y1=Atrans*y;
result.gap=y1'*y1-2*x'*y1+x'*x;
result.Dos=norm(Emax_area-Acum*y);
temp=Emax_area-Acum*y;
for i=1:NA
    result.Dos_A(i)=norm(temp((i-1)*24+1:i*24));
end
y=reshape(y,T,NA);%各分区的充电负荷
result.Cost=Price'*sum(y')';
result.Load=sum(y')+Loadprofile';
result.Loadprofile=Loadprofile;

EV_un=cumsum(EV_un);
result.EV_un=EV_un;
result.Dos_original=Emax_area-reshape(EV_un,T*NA,1);

% plot(x+Loadprofile)
% hold on;
% plot(Loadprofile,'k');
% %获得累计电量上限
% figure
% plot(std(Load_area,0,2),'k');
% hold on
% plot(std(Load_area+y,0,2));
result.y=y;
result.x=x;
result.time=toc;
%save Result