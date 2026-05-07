clear all

N=500; alpha=1; mu_cons=10; nu_cons=15;
lam=20; mu_lam=10; nu_lam=15;
N_reactants=(2:10);
samples=10;

%initialize
for i=1:length(N_reactants)
    v0{i}(1)=1;
end

parfor i=1:length(N_reactants)
    k_=repmat((1:N_reactants(i)),1,samples);
    for q=1:length(k_) %place of the sizer
        k=k_(q);
        gen=1; v0{i}(gen)=1; v=v0{i}(gen);
        x=zeros(1,N_reactants(i));
        c=x./v;
        lambda_=10.^unifrnd(log10(0.1),log10(100),1,N_reactants(i));
        %lambda_=unifrnd(1,100,1,N_reactants(i));
        mu_sample=10.^unifrnd(log10(mu_cons),log10(100));
        %mu_sample=unifrnd(mu_cons,100);
        nu_sample=10.^unifrnd(log10(nu_cons),log10(100));
        %nu_sample=unifrnd(nu_cons,100);
        c_=10*ones(1,N_reactants(i));
        rec=[ones(1,2*N_reactants(i)) -ones(1,N_reactants(i))];
        while N>gen
            t=0; v=v0{i}(gen);
            %checkpoints
            for j=1:N_reactants(i)
                mu_=zeros(1,N_reactants(i));
                nu_=zeros(1,N_reactants(i));
                if j==k
                    mu_(j)=mu_lam;
                    nu_(j)=nu_lam;
                else
                    mu_(j)=mu_sample;
                    nu_(j)=nu_sample;
                end
            while c(j)<c_(j)
                r1=rand;r2=rand;
                mu=1;
                a=[mu_.*v^2,nu_.*v,lambda_.*x];
                a0=sum(a);
                tau=(1/a0)*log(1/r1);
                t=t+tau;
                %v=v0{i}(gen)*exp(alpha*t); %exponential volume growth
                v=v0{i}(gen)+alpha*t; %linear volume growth
                while ~((sum(a(1:mu-1)) < r2*a0) && (r2*a0 <= sum(a(1:mu))))
                    mu=mu+1;
                end
                x(mod(mu-1,N_reactants(i))+1)=x(mod(mu-1,N_reactants(i))+1)+rec(mu);
                c=x./v;
            end
            end
            gen=gen+1; 
            v0{i}(gen)=v/2;
            x=floor(x/2);
            c=x./v0{i}(gen);
        end
        mdl = fitlm(v0{i}(100:end-1),v0{i}(101:end));
        slp{i}(q)=mdl.Coefficients.(1)(2);
        
        mdl = fitlm(v0{i}(100:end-1),2*v0{i}(101:end)-v0{i}(100:end-1));
        slp_{i}(q)=mdl.Coefficients.(1)(2);

        rho{i}(q)=lambda_(k)/alpha;
    end
end

col={'r','g','c','b','m',[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840]};
symb={'o', '+', '*', 'x', 'square', 'diamond', '^', '<', '>','pentagram'};
count=1;
figure(1)
clf
hold on
for i=1:length(N_reactants)
    k_=repmat((1:N_reactants(i)),1,samples);
    for q=1:length(k_) %place of the sizer
        k=k_(q);
        rho_all(count)=rho{i}(q); slp_all(count)=slp{i}(q); count=count+1;
        plot(rho{i}(q),slp{i}(q),'color',col{k},'Marker',symb{i},'linewidth',1.5,'MarkerSize',20,'HandleVisibility','off')
    end
    plot(NaN,'linestyle','none','color','k','Marker',symb{i},'linewidth',1,'MarkerSize',20,'DisplayName',"N="+N_reactants(i))
end

clear E Y rho_all_mean rho_all_std slp_all_mean slp_all_std
rho_all_logged=log10(rho_all);
[Y,E]=discretize(rho_all_logged,10);

for j=1:length(E)-1
    rho_all_mean(j)=mean(rho_all_logged(E(j)<=rho_all_logged & rho_all_logged<=E(j+1)));
    rho_all_std(j)=std(rho_all_logged(E(j)<=rho_all_logged & rho_all_logged<=E(j+1)));

    slp_all_mean(j)=mean(slp_all(E(j)<=rho_all_logged & rho_all_logged<=E(j+1)));
    slp_all_std(j)=std(slp_all(E(j)<=rho_all_logged & rho_all_logged<=E(j+1)));
end
plot(10.^rho_all_mean,slp_all_mean,'k','LineWidth',2,'HandleVisibility','off')
plot(10.^rho_all_mean,slp_all_mean+slp_all_std,'k',10.^rho_all_mean,slp_all_mean-slp_all_std,'k','HandleVisibility','off')

for q=1:max(N_reactants)
    plot(NaN,'color',col{q},'linewidth',1.5,'MarkerSize',20,'DisplayName',"k="+q)
end

ylabel("$f^{'}$",'interpreter','latex')
xlabel('$\lambda_k/\alpha$','interpreter','latex')
legend show
set(gca, 'FontSize',70,'xscale','log')
set(gca,'TickLabelInterpreter','latex','FontSize',70)
ylim([-0.25 1])
xlim([0.1 100])
pbaspect([1 1 1])
box on
