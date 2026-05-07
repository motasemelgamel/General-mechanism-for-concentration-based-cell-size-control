clear all

N=1000;
alpha=1; mu_cons=1; nu_cons=1;
lam=10; mu_lam=1; nu_lam=15;
N_reactants=3;
k=2; %place of the sizer

q=1; gen=1; v0(gen)=1; v=v0(gen); t_total(q)=0; v_total(q)=v0(gen);
x=zeros(1,N_reactants);
for j_=1:N_reactants
    c{j_}(q)=x(j_)./v;
end
lambda_=1*ones(1,N_reactants); lambda_(k)=lam;
c_=10*ones(1,N_reactants);
rec=[ones(1,2*N_reactants) -ones(1,N_reactants)];
while N>gen
    t=0; v=v0(gen);
    %checkpoints
    for j=1:N_reactants
        mu_=zeros(1,N_reactants);
        nu_=zeros(1,N_reactants);
        if j==k
            mu_(j)=mu_lam;
            nu_(j)=nu_lam;
        else
            mu_(j)=mu_cons;
            nu_(j)=nu_cons;
        end
    while c{j}(q)<c_(j)
        r1=rand;r2=rand;
        mu=1;
        a=[mu_.*v^2,nu_.*v,lambda_.*x];
        a0=sum(a);
        tau=(1/a0)*log(1/r1);
        t=t+tau; t_total(q+1)=t_total(q)+tau;
        v=v0(gen)*exp(alpha*t); v_total(q)=v;
        while ~((sum(a(1:mu-1)) < r2*a0) && (r2*a0 <= sum(a(1:mu))))
            mu=mu+1;
        end
        x(mod(mu-1,N_reactants)+1)=x(mod(mu-1,N_reactants)+1)+rec(mu);
        q=q+1;
        for j_=1:N_reactants
            c{j_}(q)=x(j_)./v;
      end
    end
    end
    T(gen)=t;
    gen=gen+1; 
    v0(gen)=v/2; v_total(q)=v/2;
    x=floor(x/2);
    for j_=1:N_reactants
        c{j_}(q)=x(j_)./v0(gen);
    end
end
mdl = fitlm(v0(200:end-1),v0(201:end));
slp=mdl.Coefficients.(1);

mdl = fitlm(v0(200:end-1),2*v0(201:end)-v0(200:end-1));
slp_=mdl.Coefficients.(1);

eps=log(v0/mean(v0));
delta=T*alpha-mean(T*alpha);
mdl = fitlm(eps(200:end-1),delta(200:end),'Intercept',false);
beta=-mdl.Coefficients.(1);

figure(1)
clf
subplot(1,2,1)
hold on
plot(v0(200:end-1)/mean(v0(200:end-1)),v0(201:end)/mean(v0(201:end)),'b.',(v0(200:end-1))/mean(v0(200:end-1)),(slp(2)*v0(200:end-1)+slp(1))/mean(v0(200:end-1)),'k','linewidth',3,'MarkerSize',35)
ylabel('$b_{n+1}/\bar{b}$','interpreter','latex')
xlabel('$b_n/\bar{b}$','interpreter','latex')
set(gca, 'FontSize',60)
xlim([0.88 1.15])
ylim([0.88 1.15])
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',60)
box on

subplot(1,2,2)
hold on
plot((0.5:0.1:1.5),0*(0.5:0.1:1.5)+1,'k','linewidth',3,'MarkerSize',25)
plot((0.5:0.1:1.5),0.5*(0.5:0.1:1.5)+0.5,'k','linewidth',3,'MarkerSize',25)
plot((0.5:0.1:1.5),1*(0.5:0.1:1.5),'k','linewidth',3,'MarkerSize',25)
ylabel('$b_{n+1}/\bar{b}$','interpreter','latex')
xlabel('$b_n/\bar{b}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',60)
set(gca, 'FontSize',60)
pbaspect([1 1 1])
box on

figure(2)
clf
subplot(1,2,1)
hold on
plot(v0(200:end-1)/mean(v0(200:end-1)),v0(201:end)/mean(v0(201:end)),'k.',(v0(200:end-1))/mean(v0(200:end-1)),(slp(2)*v0(200:end-1)+slp(1))/mean(v0(200:end-1)),'r','linewidth',3,'MarkerSize',25)
ylabel('$b_{n+1}/\bar{b}$','interpreter','latex')
xlabel('$b_n/\bar{b}$','interpreter','latex')
set(gca, 'FontSize',50)
pbaspect([1 1 1])
box on

subplot(1,2,2)
hold on
plot(v0(200:end-1)/mean(v0(200:end-1)),(2*v0(201:end)-v0(200:end-1))/mean(2*v0(201:end)-v0(200:end-1)),'k.',v0(200:end-1)/mean(v0(200:end-1)),(slp_(2)*v0(200:end-1)+slp_(1))/mean(v0(200:end-1)),'r','linewidth',3,'MarkerSize',25)
ylabel('rescaled added size, \Delta/\langle\Delta\rangle')
xlabel('rescaled birth size, $b_n/\bar{b}$','interpreter','latex')
set(gca, 'FontSize',50)
pbaspect([1 1 1])
box on

% subplot(1,2,2)
% hold on
% plot(eps(200:end-1),delta(200:end)+mean(T),'k.',eps(200:end-1),-beta*eps(200:end-1)+mean(T),'r','linewidth',3,'MarkerSize',25)
% ylabel('\alpha_n T_n')
% xlabel('ln(b_n/\langle b\rangle)')
% set(gca, 'FontSize',30)
% pbaspect([1 1 1])
% box on