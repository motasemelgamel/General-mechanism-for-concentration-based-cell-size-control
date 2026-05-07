clear all

N=500;
alpha=1; mu_cons=1; nu_cons=1;
lam=10; mu_lam=1; nu_lam=15;
N_reactants=3;
k=2; %place of the sizer

q=1; gen=1; v0(gen)=1; v=v0(gen); t_total(q)=0; v_total(q)=v0(gen);
x=zeros(1,N_reactants);
for j_=1:N_reactants
    c{j_}(q)=x(j_)./v;
    c_num{j_}(q)=x(j_)./v;
end
lambda_=1*ones(1,N_reactants); lambda_(k)=lam;
c_=10*ones(1,N_reactants);
rec=[ones(1,2*N_reactants) -ones(1,N_reactants)];
p=1;
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
        t_mol=0; b_mol{j}(gen)=v;
    while c{j}(q)<c_(j)
        r1=rand;r2=rand;
        mu=1;
        a=[mu_.*v^2,nu_.*v,lambda_.*x];
        a0=sum(a);
        tau=(1/a0)*log(1/r1);
        t=t+tau; t_total(q+1)=t_total(q)+tau;
        t_mol=t_mol+tau;
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
    T_mol{j}(gen)=t_mol;
    if t_mol<=10/lambda_(k) && j==k
        small_gen(p)=gen;
        p=p+1;
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

for j=1:length(T_mol)
    t_mu(j)=mean(T_mol{j});
end

%numerical
dt=0.0001; gen=1; time=0;
q=1; v0_num(gen)=1; v_num=v0_num(gen);

while N>gen
    t=0;
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
        %c0_(1)=c_num{1}(end); c0_(2)=c_num{2}(end);
        while c_num{j}(q)<c_(j)
            for p=1:N_reactants
                c_num{p}(q+1)=c_num{p}(q) + (mu_(p)*v_num+nu_(p)-(lambda_(p)+alpha)*c_num{p}(q))*dt;
                %c_th{p}(q)=c0_(p).*exp(-(alpha+lambda_(p)).*t) + exp(-(alpha+lambda_(p)).*t).*((v0_num(end)*mu_(p)/(2*alpha+lambda_(p))).*(exp((2*alpha+lambda_(p)).*t)-1) +(nu_(p)./(alpha+lambda_(p))).*(exp((alpha+lambda_(p)).*t)-1) );
            end
            time=time+dt;
            t=t+dt;
            v_num=v0_num(gen)*exp(alpha*t);
            q=q+1;
        end
    end
    gen=gen+1;
    v0_num(gen)=v_num/2;
    v_num=v0_num(gen);
end


lim=[sum(T(1:300)) sum(T(1:306))];
figure(1)
clf
sgt=sgtitle("Number of stages, N="+N_reactants+", sizer stage, k="+k,'interpreter','latex');
sgt.FontSize = 40;
subplot(2,6,1:5)
hold on
plot(t_total,c{1},'b',[0 max(t_total)],c_(1)*[1 1],'b--',t_total,c{2},'r-',[0 max(t_total)],c_(2)*[1 1],'r--',t_total,c{3},'g',[0 max(t_total)],c_(3)*[1 1],'k--','linewidth',2,'MarkerSize',25)
plot((0:dt:time),c_num{1}(1:end-1),'b-.',(0:dt:time),c_num{2}(1:end-1),'r-.',(0:dt:time),c_num{3}(1:end-1),'g-.','linewidth',2,'MarkerSize',15)
ylabel('Concentration','interpreter','latex')
set(gca, 'FontSize',40)
h(1)=plot(NaN,NaN,'k-.','MarkerSize',10,'linewidth',1.5);
h(2)=plot(NaN,NaN,'k','MarkerSize',10,'linewidth',1.5);
h(3)=plot(NaN,NaN,'b','MarkerSize',10,'linewidth',1.5);
h(4)=plot(NaN,NaN,'r','MarkerSize',10,'linewidth',1.5);
h(5)=plot(NaN,NaN,'g','MarkerSize',10,'linewidth',1.5);
legend(h,'numerical','stochastic','k=1','k=2','k=3','interpreter','latex')
xlim(lim)
ylim([0 11])
xticks([sum(T(1:300)) sum(T(1:300))+1 sum(T(1:300))+2 sum(T(1:300))+3 sum(T(1:300))+4])
yticks([0 10])
yticklabels({'0','$c^*$'})
xticklabels(0:1:4)
set(gca,'TickLabelInterpreter','latex','FontSize',50)

subplot(2,6,7:11)
plot(t_total,v_total/(mean(v0)),'k','linewidth',2,'MarkerSize',25)
ylabel('Size, S/$\bar{b}$','interpreter','latex')
xlabel('Time, $\alpha t$','interpreter','latex')
set(gca, 'FontSize',40)
xlim(lim)
xticks([sum(T(1:300)) sum(T(1:300))+1 sum(T(1:300))+2 sum(T(1:300))+3 sum(T(1:300))+4])
xticklabels(0:1:4)
set(gca,'TickLabelInterpreter','latex','FontSize',50)