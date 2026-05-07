clear all

N=500;
alpha=1; mu_cons=1; nu_cons=1;
N_reactants=2;

mu_arr=logspace(log10(1e-3),log10(1),30);
c_th=[0.5 0.8 1];
q_=1; %altered sizer or non-sizer
for m=1:length(c_th)
for k=1:N_reactants %place of the sizer
for i=1:length(mu_arr)
    clear v0 c
    gen=1; v0(gen)=1; v=v0(gen);
    x=zeros(1,N_reactants);
    for j_=1:N_reactants
        c{j_}=x(j_)./v;
    end
    lambda_=ones(1,N_reactants); lambda_(k)=10;
    c_=c_th(m)*ones(1,N_reactants); c_(k)=10;
    rec=[ones(1,2*N_reactants) -ones(1,N_reactants)];

while N>gen
    t=0; v=v0(gen);
    %checkpoints
    for j=1:N_reactants
        mu_=zeros(1,N_reactants);
        nu_=zeros(1,N_reactants); nu_(j)=nu_cons;
        if q_==1 %alter non-sizer
            if j==k 
                mu_(j)=mu_cons;
            else
                mu_(j)=mu_arr(i);
            end
        elseif q_==2 %alter sizer
            if j~=k
                mu_(j)=mu_cons;
            else
                mu_(j)=mu_arr(i);
            end
        end
    while c{j}<c_(j)
        r1=rand;r2=rand;
        mu=1;
        a=[mu_.*v^2,nu_.*v,lambda_.*x];
        a0=sum(a);
        tau=(1/a0)*log(1/r1);
        t=t+tau;
        v=v0(gen)*exp(alpha*t);
        while ~((sum(a(1:mu-1)) < r2*a0) && (r2*a0 <= sum(a(1:mu))))
            mu=mu+1;
        end
        x(mod(mu-1,N_reactants)+1)=x(mod(mu-1,N_reactants)+1)+rec(mu);
        for j_=1:N_reactants
            c{j_}=x(j_)./v;
        end
    end
    end
    gen=gen+1; 
    v0(gen)=v/2;
    x=floor(x/2);
    for j_=1:N_reactants
        c{j_}=x(j_)./v0(gen);
    end
end

mdl = fitlm(v0(100:end-1),v0(101:end));
slp=mdl.Coefficients.(1);
f{k,m}(i)=slp(2);
b_bar{k,m}(i)=mean(v0(100:end-1));
cv{k,m}(i)=std(v0(100:end-1))./mean(v0(100:end-1));
% f{k,q_}(i)=slp(2);
% b_bar{k,q_}(i)=mean(v0(100:end-1));
% cv{k,q_}(i)=std(v0(100:end-1))./mean(v0(100:end-1));
end
end
end

figure(1)
clf
subplot(1,3,1)
hold on
plot(1./mu_arr,f{1,1},'b.',1./mu_arr,f{2,1},'r.','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,f{1,2},'b^',1./mu_arr,f{2,2},'r^','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,f{1,3},'b*',1./mu_arr,f{2,3},'r*','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,f{1,4},'b+',1./mu_arr,f{2,4},'r+','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,f{1,2},'k.',1./mu_arr,f{2,2},'b.','MarkerSize',70,'LineWidth',2)
%plot(1./mu_arr,f{1,1},'g.',1./mu_arr,f{2,1},'r.','MarkerSize',70,'LineWidth',2)
plot(NaN,'color','k','Marker','.','linewidth',2,'MarkerSize',50,'DisplayName',"c*="+c_th(1))
plot(NaN,'color','k','Marker','^','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(2))
plot(NaN,'color','k','Marker','*','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(3))
%plot(NaN,'color','k','Marker','+','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(4))
set(gca,'FontSize',48,'xscale','log')
xlabel('$1/\mu$','interpreter','latex')
ylabel('$f^{''}$','interpreter','latex')
%legend('k=1, vary \mu','k=2, vary \mu')
%legend('k=1, vary \mu','k=2, vary \mu','k=1, vary \mu_s','k=2, vary \mu_s')
legend show
pbaspect([1 1 1])
ylim([-0.4 1])
box on

subplot(1,3,2)
hold on
plot(1./mu_arr,b_bar{1,1}/min(b_bar{1,1}),'b.',1./mu_arr,b_bar{2,1}/min(b_bar{2,1}),'r.','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,b_bar{1,2}/min(b_bar{1,2}),'b^',1./mu_arr,b_bar{2,2}/min(b_bar{2,2}),'r^','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,b_bar{1,3}/min(b_bar{1,3}),'b*',1./mu_arr,b_bar{2,3}/min(b_bar{2,3}),'r*','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,b_bar{1,4}/min(b_bar{1,4}),'b+',1./mu_arr,b_bar{2,4}/min(b_bar{2,4}),'r+','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,b_bar{1,2}/min(b_bar{1,2}),'k.',1./mu_arr,b_bar{2,2}/min(b_bar{2,2}),'b.','MarkerSize',70,'LineWidth',2)
%plot(1./mu_arr,b_bar{1,1},'g.',1./mu_arr,b_bar{2,1},'r.','MarkerSize',70,'LineWidth',2)
plot(NaN,'color','k','Marker','.','linewidth',2,'MarkerSize',50,'DisplayName',"c*="+c_th(1))
plot(NaN,'color','k','Marker','^','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(2))
plot(NaN,'color','k','Marker','*','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(3))
%plot(NaN,'color','k','Marker','+','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(4))
set(gca,'FontSize',48,'xscale','log')
xlabel('$1/\mu$','interpreter','latex')
ylabel('$\bar{b}/b_{min}$','interpreter','latex')
%legend('k=1, vary \mu','k=2, vary \mu')
%legend('k=1, vary \mu','k=2, vary \mu','k=1, vary \mu_s','k=2, vary \mu_s')
legend show
pbaspect([1 1 1])
ylim([0.9 3])
box on

subplot(1,3,3)
hold on
plot(1./mu_arr,cv{1,1},'b.',1./mu_arr,cv{2,1},'r.','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,cv{1,2},'b^',1./mu_arr,cv{2,2},'r^','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
plot(1./mu_arr,cv{1,3},'b*',1./mu_arr,cv{2,3},'r*','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,cv{1,4},'b+',1./mu_arr,cv{2,4},'r+','MarkerSize',15,'LineWidth',2,'HandleVisibility','off')
%plot(1./mu_arr,cv{1,2},'k.',1./mu_arr,cv{2,2},'b.','MarkerSize',70,'LineWidth',2)
%plot(1./mu_arr,cv{1,1},'g.',1./mu_arr,cv{2,1},'r.','MarkerSize',70,'LineWidth',2)
plot(NaN,'color','k','Marker','.','linewidth',2,'MarkerSize',50,'DisplayName',"c*="+c_th(1))
plot(NaN,'color','k','Marker','^','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(2))
plot(NaN,'color','k','Marker','*','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(3))
%plot(NaN,'color','k','Marker','+','linewidth',2,'MarkerSize',15,'DisplayName',"c*="+c_th(4))
set(gca,'FontSize',48,'xscale','log')
xlabel('$1/\mu$','interpreter','latex')
ylabel('$CV$','interpreter','latex')
%legend('k=1, vary \mu','k=2, vary \mu')
%legend('k=1, vary \mu','k=2, vary \mu','k=1, vary \mu_s','k=2, vary \mu_s')
legend show
pbaspect([1 1 1])
ylim([0 0.26])
box on

%Data
b_bar_exp{1,1}=12.4; cv_exp{1,1}=0.09;
b_bar_exp{1,2}=16.8; cv_exp{1,2}=0.25;
b_bar_exp{2,1}=19; cv_exp{2,1}=0.15;

figure(2)
clf
subplot(1,2,1)
hold on
plot([0 1],[b_bar_exp{1,1} b_bar_exp{1,2}]./b_bar_exp{1,1},'k-',[0 1],[b_bar_exp{1,1} b_bar_exp{2,1}]./b_bar_exp{1,1},'k-','MarkerSize',40,'LineWidth',2)
plot([0 0],[b_bar_exp{1,1} b_bar_exp{1,1}]./b_bar_exp{1,1},'k.',[1 1],[b_bar_exp{1,2} b_bar_exp{1,2}]./b_bar_exp{1,1},'ks',[1 1],[b_bar_exp{2,1} b_bar_exp{2,1}]./b_bar_exp{1,1},'k^','MarkerSize',40,'LineWidth',2)
set(gca,'FontSize',70,'xscale','linear')
ylabel('$\bar{b}/b_{min}$','interpreter','latex')
xlim([-0.1 1.1])
xticks([0 1])
xticklabels({'WT','\Delta'})
legend('﻿yFS110 wt','','','yFS1123 cdc13-ubi-\Delta51-70','yFS887 nmt41:cdc25')
box on
pbaspect([1 1 1])
ylim([0.9 1.6])

subplot(1,2,2)
hold on
plot([0 1],[cv_exp{1,1} cv_exp{1,2}],'k-',[0 1],[cv_exp{1,1} cv_exp{2,1}],'k-','MarkerSize',40,'LineWidth',2)
plot([0 0],[cv_exp{1,1} cv_exp{1,1}],'k.',[1 1],[cv_exp{1,2} cv_exp{1,2}],'ks',[1 1],[cv_exp{2,1} cv_exp{2,1}],'k^','MarkerSize',40,'LineWidth',2)
set(gca,'FontSize',70,'xscale','linear')
ylabel('$CV$','interpreter','latex')
xlim([-0.1 1.1])
xticks([0 1])
xticklabels({'WT','\Delta'})
box on
pbaspect([1 1 1])
ylim([0.05 0.26])