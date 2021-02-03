%% params_input
params.c=2.99792458e8;
params.e=1.60218e-19;
params.mp=1.6726231e-27;
params.me=9.1093897e-31;
params.m0=params.me/params.mp;
params.bz = 5.4;  
params.n0e = 57.5;
params.k=0.3125;
params.delta = 0.8; 
params.omega=0;
params.lambda=0.83;
params.te_=38;
params.lambda1=180/params.te_;
params.te2=params.te_*params.lambda1*params.e*1e7;


params.k1=0;
params.k2=2.5e5;
params.k3=0;
params.k4=0;
params.tau_=310/params.te_;
params.tau=310/params.te_/params.lambda;
params.te = params.te_*params.lambda;
params.ti_ = params.te_*params.tau_;
params.ti = params.te*params.tau;
params.z=1;
params.m = params.m0/params.z;
params.T0 = (params.te_*params.lambda*(params.z+params.tau))*params.e*1e7;
params.B0 = 1e5*sqrt(4*pi*params.n0e*params.T0); 
params.be=2.9/params.B0;
params.be_=-1.1/params.B0;
params.n0i = (params.n0e*(1-params.delta)/(1-params.be*params.B0/params.bz)*exp(params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7) +...
    params.n0e*params.delta/(1-params.be_*params.B0/params.bz));

params.di = sqrt(2.998e10^2*params.mp*1e3/(4*pi*params.n0e*(2.998e9*1.6e-19)^2));
params.D=params.di*params.k;
params.Omega0=-1/sqrt((params.D)^4/(3*10^10)^2/(params.T0)*4*pi*(params.e*3*10^9)^2*params.n0e);
params.Omegai = params.Omega0 * params.omega/(params.omega-1);
params.Omegae = params.Omega0/(params.omega-1);
params.ae=params.m/(params.m + params.omega^2);
params.ai=params.omega^2/(params.m + params.omega^2);
params.epsilon = (params.omega^2+params.m)*params.z/params.k^2/(params.omega-1)^2;

params.te2=params.te_*params.lambda1*params.e*1e7;
params.nita=params.c*1e2*params.T0/params.e/params.c/10/params.Omega0;
params.v0 = sqrt(params.te*params.e*2/params.me);
params.v0i = sqrt(params.ti*params.e*2/params.mp);
params.qujian=1.5;
%% load 初始电磁场
inputs=[0:0.01:50];
sol = ode45(@(x,y)vpt2(x,y,params),inputs, [0;-params.bz/params.B0;0;0;0]);
[y, yp] = deval(inputs,sol);
x = sqrt(2*inputs);
E = -yp(5,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100;
B = -y(2,:)*params.B0;
Bphi=-params.nita*yp(3,:).*x/params.D*1e5;
data={y,x,E,B,Bphi};
disp('initial state finished')
%% ode+电磁场存储
jieshumax=4;
for i=1:jieshumax
    a=max(data{i,2});
    a=a-0.9;
    a=a^2/2;
    inputs=[0:0.001:0.1,0.11:0.01:0.2,0.25:0.05:a];
    sol = ode45(@(xx,y)vptl3(xx,y,data{i,2},data{i,3},data{i,4},params),inputs, [0;-params.bz/params.B0;0;0;0]);
    [y, yp] = deval(inputs,sol);
    x = sqrt(2*inputs);
    E = -yp(5,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100;
    B = -y(2,:)*params.B0;
    Bphi=-params.nita*yp(3,:).*x/params.D*1e5;
    data{i+1,1}=y;
    data{i+1,2}=x;
    data{i+1,3}=E;
    data{i+1,4}=B;
    data{i+1,5}=Bphi;
    disp(i)
end
disp('ode finished')
data{6,1}=params;
save data_end data
%%
jieshu=jieshumax;
x=data{jieshu,2};
E=data{jieshu,3};
B=data{jieshu,4};
y=data{jieshu+1,1};
xnew=data{jieshu+1,2};
Enew=data{jieshu+1,3};
Bnew=data{jieshu+1,4};
qujian=params.qujian;
v0=params.v0;
v0i=params.v0i;
params.v0i = sqrt(params.ti*params.e*2/params.mp);
rr=@(xx,yy)sqrt(xx.^2+yy.^2);
E_e=@(xx)interp1(x,E,xx,'spline');
B_e=@(xx)interp1(x,B,xx,'spline');
        rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
        rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
        rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
        r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calcute miu
        vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
        vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
        B_tidu=dif2(B,x)/params.D*100*1e-9;
        B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
        vD=@(vx,vy,xx,yy)1/2*params.me*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
            .*B_tidu_e(r_(vx,vy,xx,yy))*1e18;
        vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
        vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
miu_guiyi1=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
miu_guiyi2=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be_*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
f_e_rec_1= @(vx,vy,xx,yy)params.n0e*(1-params.delta)*(1/pi)*exp(...
        -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy) ...
        + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
        (1+params.tau)*(interp1(xnew,y(5,:),rr(xx,yy),'spline')-params.k2*1e2/params.Omega0*(interp1(xnew,y(3,:),rr(xx,yy),'spline')) ...
        -1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline'))+params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7); 
f_e_rec_2=@(vx,vy,xx,yy)params.n0e*params.delta*(params.lambda/pi)*exp(...
        -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy)) ...
        + params.lambda*(1+params.tau)*interp1(xnew,y(5,:),rr(xx,yy),'spline'));
f_e_rec__=@(vx,vy,xx,yy)f_e_rec_1(vx,vy,xx,yy)+f_e_rec_2(vx,vy,xx,yy);
nn_e=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-qujian,qujian,-qujian,qujian);
params.n0i=nn_e(0,0);
f_i_rec__=@(vx,vy,xx,yy)params.n0i*(1-params.delta)*(1/pi)^1*exp(...
        -(vx.^2+vy.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(5,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^1*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(5,:),rr(xx,yy),'spline')); 
    
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-10,10,-10,10); 
%% function    
function dydx=vptl3(xx,y,x,E,B,params)
        dydx=zeros(5,1);  
        lambda = params.lambda;
        tau = params.tau;
        tau_ = params.tau_;
        delta = params.delta;
        k = params.k;
        m = params.m/1;
        z = params.z;
        be = params.be;
        be_ = params.be_;
        qujian=params.qujian;
        v0 = sqrt(params.te*1.6e-19*2/params.me);
        rr=@(xx,yy)sqrt(xx.^2+yy.^2);
        E_e=@(xx)interp1(x,E,xx,'spline');
        B_e=@(xx)interp1(x,B,xx,'spline');
        rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
        rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
        rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
        r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calcute miu
        vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
        vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
        B_tidu=dif2(B,x)/params.D*100*1e-9;
        B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
        vD=@(vx,vy,xx,yy)1/2*params.me.*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
            .*B_tidu_e(r_(vx,vy,xx,yy))*1e18; %./sqrt(1-(vx.^2+vy.^2).*v0^2./params.c.^2)
        vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
        vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        miu_guiyi1=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be./(-y2)/v0^2;
        miu_guiyi2=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be_./(-y2)/v0^2;
          
        f_e_rec_1= @(vx,vy,xx,yy,y2,y1,y5,y3)params.n0e*(1-params.delta)*(1/pi)*exp(...
            -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy,y2) ...
            + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
            (1+params.tau)*(y5-params.k2*1e2/params.Omega0*y3-1/(params.omega-1)*y1)+params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7); 
        f_e_rec_2= @(vx,vy,xx,yy,y2,y1,y5)params.n0e*params.delta*(params.lambda/pi)*exp(...
            -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy,y2)) ...
            + params.lambda*(1+params.tau)*y5);
        f_e_rec__=@(vx,vy,xx,yy,y2,y1,y5,y3)f_e_rec_1(vx,vy,xx,yy,y2,y1,y5,y3)+f_e_rec_2(vx,vy,xx,yy,y2,y1,y5);
    Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    params.n0i=integral2(@(vx,vy)f_e_rec__(vx,vy,1e-7,0,-params.bz/params.B0,0,0,0),-qujian,qujian,-qujian,qujian);
    nie = params.n0i/params.n0e;
    Ee_func=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*t),0,y(2),y(1),y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y2=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,t,y(1),y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y1=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),t,y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y5=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),t,y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y3=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),t),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_=integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*t),0,y(2),y(1),y(5)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y5=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),t),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y2=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,t,y(1),y(5)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*xx+(params.omega/(params.omega-1))*y(1)-y(5)));
    Ei_=exp(-lambda*(1+tau)/tau_*y(5));
    dydx(1)=y(2);
    if xx<1e-5
        xx_=1e-5;
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx_),0,y(2),y(1),y(5),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
        Ee_=integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx_*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        dydx(3)=y(4);
        nita=params.nita;
        dydx(4)=(params.D^2/nita*4*pi/params.c*params.k2*params.n0e*(1-delta)*Ee*params.e*params.c*10-2*y(4))/2/(xx+1e-7);
        
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7;
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta);
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);  
        Phi=-(Ee_func_y3(y(3)+1e-7)-Ee)/1e-7;
        K=(Ee_func_y5(y(5)+1e-7)-Ee)/1e-7+(Ee__func_y5(y(5)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_;  
        dydx(5) = (Psi*y(2)+PPsi*dydx(2)+Theta+Phi*y(4))/K;
    else
        
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
        Ee_=integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        while abs(dydx(2))>1
            je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx),0,y(2),y(1),y(5),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
            dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
            qujian=qujian-0.1;
        end
        dydx(3)=y(4);
        nita=params.nita;
        dydx(4)=(params.D^2/nita*4*pi/params.c*params.k2*params.n0e*(1-delta)*Ee*params.e*params.c*10-2*y(4))/2/(xx+1e-7);
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7;
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta);
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);  
        Phi=-(Ee_func_y3(y(3)+1e-10)-Ee)/1e-10;
        K=(Ee_func_y5(y(5)+1e-7)-Ee)/1e-7+(Ee__func_y5(y(5)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_;  
        dydx(5) = (Psi*y(2)+PPsi*dydx(2)+Theta+Phi*y(4))/K;
    end    
end
function dydx=vpt2(x,y,params)
    % params: 常数结构体，含有下列参数
    lambda = params.lambda;
    tau = params.tau;
    tau_ = params.tau_;
    delta = params.delta;
    
    k = params.k;
    m = params.m/1;
    z = params.z;
    be = params.be;
    be_ = params.be_;
    nita=params.nita;
    nie = params.n0i/params.n0e;
    
    dydx=zeros(5,1);
    xie=be/y(2);
    xie_ = be_/y(2);
    Ee=exp((z+tau)*(params.ae*params.epsilon*x/(1+xie)- ...
        (1/(params.omega-1))*y(1)+y(5)-params.k2*1e2/params.Omega0*y(3)) ...
        +params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7);
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*x+(params.omega/(params.omega-1))*y(1)-y(5)));
   
    Ee_=exp(lambda*(1+tau)*y(5));
    Ei_=exp(-lambda*(1+tau)/tau_*y(5));
    
    dydx(1)=y(2);
    dydx(2)=(1-delta)*(Ee/(params.omega-1)/(1+xie)^2 - nie*Ei*params.omega/(params.omega-1));
    dydx(3)=y(4);
    dydx(4)=(params.D^2/nita*4*pi/params.c*params.k2*params.n0e*(1-delta)/(1+xie)*Ee*params.e*params.c*10-2*y(4))/2/(x+1e-7);
   
    
    Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei+Ee/(1+xie));
    PPsi=-Ee*(1+(1+tau)*params.ae*params.epsilon*x/(1+xie))*(be/(be+y(2))^2)...
        -Ee_*delta/(1-delta)*(be_/(be_+y(2))^2);
    Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau-params.ae*Ee/(1+xie)^2);
    Phi=(1+tau)/(1+xie)*params.k2*1e2/params.Omega0*Ee;
    K=(1+tau)*(Ee/(1+xie)+lambda*delta/(1-delta)/(1+xie_)*Ee_+nie/tau*Ei+nie*lambda/tau_*delta/(1-delta)*Ei_);
    dydx(5) = (Psi*y(2)+PPsi*dydx(2)+Theta+Phi*y(4))/K;
end

function result=judge_0(x,y)
    if x==0 & y==0
        result=0;
    else
        result=1./sqrt(x.^2+y.^2);
    end
end