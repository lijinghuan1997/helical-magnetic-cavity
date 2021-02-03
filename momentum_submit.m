%% 
% load data
v0 = params.v0;
v0i = params.v0i;
jieshumax=4;
jieshu=jieshumax;
x=data{jieshu,2};
E=data{jieshu,3};
B=data{jieshu,4};

y=data{jieshu+1,1};
xnew=data{jieshu+1,2};
Enew=data{jieshu+1,3};
Bnew=data{jieshu+1,4};
Bphi=data{jieshu+1,5};
qujian=params.qujian; % integral-span
% handles
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
% distribution functions 2-D in velocity space 
f_e_rec_1= @(vx,vy,xx,yy)params.n0e*(1-params.delta)*(1/pi)*exp(...
            -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy) ...
            + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
       (1+params.tau)*(interp1(xnew,y(5,:),rr(xx,yy),'spline')-params.k2*1e2/params.Omega0*interp1(xnew,y(3,:),rr(xx,yy),'spline') ...
            -1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')) ...
            +params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7); 
f_e_rec_2= @(vx,vy,xx,yy)params.n0e*params.delta*(params.lambda/pi)*exp(...
            -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy)) ...
            + params.lambda*(1+params.tau)*interp1(xnew,y(5,:),rr(xx,yy),'spline'));
f_e_rec__=@(vx,vy,xx,yy)f_e_rec_1(vx,vy,xx,yy)+f_e_rec_2(vx,vy,xx,yy);
nn_e_1=@(xx,yy)integral2(@(vx,vy)f_e_rec_1(vx,vy,xx,yy),-qujian,qujian,-qujian,qujian);
nn_e_2=@(xx,yy)integral2(@(vx,vy)f_e_rec_2(vx,vy,xx,yy),-qujian,qujian,-qujian,qujian);
nn_e=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-qujian,qujian,-qujian,qujian);
params.n0i=nn_e(1e-7,0);    


f_i_rec__=@(vx,vy,xx,yy)params.n0i*(1-params.delta)*(1/pi)^1*exp(...
        -(vx.^2+vy.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(5,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^1*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(5,:),rr(xx,yy),'spline')); 
   
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-10,10,-10,10);
nn_=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-qujian,qujian,-qujian,qujian);
vvreal= @(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,rr(xx,yy),0).*vy*v0,-qujian,qujian,-qujian,qujian,'AbsTol',1e2,'RelTol',1e-4)/nn_e(rr(xx,yy),0);
vvreali= @(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,rr(xx,yy),0).*vy*v0i,-qujian,qujian,-qujian,qujian)/nn_i(rr(xx,yy),0);

%% calculate
xplot=[0:0.05:3.5];
%xplot=[0:0.01:3.5];
for i =1:length(xplot)
    Ne1(i)=nn_e_1(xplot(i),0);
    Ne2(i)=nn_e_2(xplot(i),0);
    Ne(i)=Ne1(i)+Ne2(i);
    Ni(i)=nn_i(xplot(i),0);
    int_ve(i)=vvreal(xplot(i),0);
    int_vi(i)=vvreali(xplot(i),0);
    Pe(i)=integral2(@(vx,vy)f_e_rec__(vx,vy,xplot(i),0).*vx.^2*v0^2,-5,5,-5,5,'AbsTol',1e10,'RelTol',1e-4)*params.me*10^6;
    Teperp(i)=Pe(i)/integral2(@(vx,vy)f_e_rec__(vx,vy,xplot(i),0),-5,5,-5,5)/1e6/params.e;
    Tepara(i)=(params.te_*params.lambda*Ne1(i)+params.te_*Ne2(i))/Ne(i);
end  
int_vz=params.k2*Ne1./Ne/1e3;

%% example plot
xwidth=0.7;
ywidth=0.1;
% magnetic field (nT)
subplot('position',[0.15 0.97-ywidth xwidth ywidth]);
plot(xnew*params.D/1e5,Bnew)
ylim([0,20])
set(gca,'xtick',[])
set(gca,'ytick',[5 15])
ylabel('B/nT')
% electric field (V/m)
subplot('position',[0.15 0.97-2*ywidth xwidth ywidth]);
plot(xnew*params.D/1e5,Enew)
ylim([-1.5e-3, 0.5e-3])
set(gca,'xtick',[])
set(gca,'ytick',[-1e-3 0 1e-3])
ylabel('E/ V/m')
% number density cm^-3
subplot('position',[0.15 0.97-3*ywidth xwidth ywidth]);
plot(xplot*params.D/1e5,Ne)
hold on
plot(xplot*params.D/1e5,Ni)
ylim([50.5,53.5])
set(gca,'xtick',[])
set(gca,'ytick',[51.5 52.5])
ylabel('N/cm^-3')
% electron temperatures eV
subplot('position',[0.15 0.97-4*ywidth xwidth ywidth]);
plot(xplot*params.D/1e5,Teperp)
hold on
plot(xplot*params.D/1e5,Tepara)
ylim([32,42])
set(gca,'xtick',[])
set(gca,'ytick',[35 39])
ylabel('Te/eV')
% electron bulk velcity km/s
subplot('position',[0.15 0.97-5*ywidth xwidth ywidth]);
plot(xplot*params.D/1e5,int_ve/1e3)
hold on
plot(xplot*params.D/1e5,int_vz)
ylim([-0.5e2, 1.5e2])
set(gca,'xtick',[])
set(gca,'ytick',[0 1e2])
ylabel('Ve/ km/s')
% Bphi  nT
subplot('position',[0.15 0.97-6*ywidth xwidth ywidth]);
plot(xnew*params.D/1e5,Bphi)
ylim([-3,1])
set(gca,'xtick',[])
set(gca,'ytick',[-2 0])
ylabel('Bphi/nT')
%% 3-D velocity distribution
f_e_1= @(vx,vy,vz,xx,yy)params.n0e*(1-params.delta)*(1/pi)^1.5*exp(...
            -vx.^2 - vy.^2  - (vz-params.k2/v0).^2 + miu_guiyi1(vx,vy,xx,yy) ...
            + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
            (1+params.tau)*(interp1(xnew,y(5,:),rr(xx,yy),'spline')-params.k2*1e2/params.Omega0*interp1(xnew,y(3,:),rr(xx,yy),'spline') ...
            -1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')) ...
            +params.me*1e3/2*params.k2^2*1e4/params.te/params.e/1e7); 
f_e_2= @(vx,vy,vz,xx,yy)params.n0e*params.delta*(params.lambda/pi)^1.5*exp(...
            -params.lambda*(vx.^2+vy.^2+vz.^2-miu_guiyi2(vx,vy,xx,yy)) ...
            + params.lambda*(1+params.tau)*interp1(xnew,y(5,:),rr(xx,yy),'spline'));
f_e_=@(vx,vy,vz,xx,yy)f_e_1(vx,vy,vz,xx,yy)+f_e_2(vx,vy,vz,xx,yy);

f_e=@(vx,vy,vz,xx)f_e_(vx,vy,vz,xx,0);
%% functions 
function result=judge_0(x,y)
    if x==0 & y==0
        result=0;
    else
        result=1./sqrt(x.^2+y.^2);
    end
end
