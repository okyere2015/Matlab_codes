function [t,S,E,I,C,D,R,S1,E1,I1,C1,D1,R1,S4,E4,I4,C4,D4,R4]=MONKEYPOX_MODEL(neta,omega,gamma,beta,mu,delta,alpha,rho,N,u,u3,S0,E0,I0,C0,D0,R0,S10,E10,I10,C10,D10,R10,S40,E40,I40,C40,D40,R40,MaxTime)

% It is the Covid-19 model.
if nargin==0
    neta=0.016744;
    omega=29.08;
    beta=0.022325;
    mu=0.00004252912;
    rho=0.036246;
    alpha=0.5;
    delta=0.003286;
    gamma=0.088366;
    
    u=0.5;
    u3=0.3;
    
    
     
   S0=30799995;
    E0=0;
    I0=5;
    C0=0;
    D0=0.0;
    R0=0.0;
    
    S10=30799995;
    E10=0;
    I10=5;
    C10=0;
    D10=0.0;
    R10=0.0;
    
     S40=30799995;
    E40=0;
    I40=5;
    C40=0;
    D40=0.0;
    R40=0.0;
    
    
    R01=(((beta*neta)./((gamma+alpha+mu+delta)*(neta+mu))))^(0.5)
   B=((0.5)*((beta*neta)./((gamma+alpha+mu+delta)*(neta+mu)))^(-0.5)).*((neta)./((gamma+alpha+mu+delta)*(neta+mu)));
   B1=((0.5)*((beta*neta)./((gamma+alpha+mu+delta)*(neta+mu)))^(-0.5))*(((neta+mu)*(gamma+delta+mu+alpha)*beta-neta*beta*(gamma+alpha+mu+delta))./((gamma+alpha+mu+delta)*(neta+mu))^2);
    B2=((0.5)*((beta*neta)./((gamma+alpha+mu+delta)*(neta+mu)))^(-0.5))*(((-beta*neta)*(neta+mu))./((gamma+alpha+mu+delta)*(neta+mu))^2);
    B3=((0.5)*((beta*neta)./((gamma+alpha+mu+delta)*(neta+mu)))^(-0.5))*(((-beta*neta)*(gamma+alpha+delta+2*mu))./((gamma+alpha+mu+delta)*(neta+mu))^2); 
   T1=B*(beta./R01)
   T2=B1*(neta./R01)
   T3=B2*(gamma./R01)
   T4=B2*(alpha./R01)
   T5=B2*(delta./R01)
   T6=B3*(mu./R01) 
    MaxTime=120;
    
  
end

S=S0;E=E0;I=I0;C=C0;D=D0;R=R0;S1=S10;E1=E10;I1=I10;C1=C10;D1=D10;R1=R10;
S4=S40;E4=E40;I4=I40;C4=C40;D4=D40;R4=R40;
N=S+E+I+C+D+R;
% The main iteration
options=odeset('RelTol',1e-3);
[t,pop]=ode45(@Diff_2_6,[0 MaxTime],[S E I C D R S1 E1 I1 C1 D1 R1 S4 E4 I4 C4 D4 R4],options,[neta omega gamma beta mu delta rho alpha N u u3]);

S=pop(:,1);E=pop(:,2);I=pop(:,3);C=pop(:,4);D=pop(:,5);R=pop(:,6);S1=pop(:,7);E1=pop(:,8);I1=pop(:,9);C1=pop(:,10);D1=pop(:,11);R1=pop(:,12);
S4=pop(:,13);E4=pop(:,14);I4=pop(:,15);C4=pop(:,16);D4=pop(:,17);R4=pop(:,18);
 % plots the graphs with scaled colours
 figure(1)
r=plot(t,S,'-b',t,S1,'.m');
legend(r,'u_1=0,u_2=0','u_1 \neq 0,u_2=0')
xlabel 'Time (days)'
ylabel 'H_S(t)'

figure(2)
V=plot(t,E,'-b',t,E1,'.m');
legend(V,'u_1=0,u_2=0','u_1 \neq 0,u_2=0')
xlabel 'Time (days)'
ylabel 'H_E (t)'

figure(3)
k=plot(t,I,'-b',t,I1,'.m');
legend(k,'u_1=0,u_2=0','u_1 \neq 0,u_2=0')
xlabel 'Time (days)'
ylabel 'H_I(t)'

figure(4)
c=plot(t,C,'-c',t,C1,'-b');
legend(c,'Hospitalised')
xlabel 'Time (days)'
ylabel 'H_Q(t)'

figure(6)
Z=plot(t,R,'-b',t,R1,'.m');
legend(Z,'u_1=0,u_2=0','u_1 \neq 0,u_2=0')
xlabel 'Time (days)'
ylabel 'H_R(t)'

figure(7)
z=plot(t,D,'-b',t,D1,'.m');
legend(z,'u_1=0,u_2=0','u_1 \neq 0,u_2=0')
xlabel 'Time (days)'
ylabel 'D(t)'

 figure(8)
L=plot(t,S,'-b',t,S4,'-c');
legend(L,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_S(t)'

figure(9)
l=plot(t,E,'-b',t,E4,'-c');
legend(l,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_E (t)'

figure(10)
K=plot(t,I,'-b',t,I4,'.c');
legend(K,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_I(t)'

figure(11)
w=plot(t,C,'-b',t,C4,'.c');
legend(w,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_Q(t)'



figure(12)
W=plot(t,R,'-b',t,R4,'-c');
legend(W,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_R(t)'

figure(13)
P=plot(t,D,'-b',t,D4,'-c');
legend(P,'u_1=0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'D(t)'

 figure(14)
Y=plot(t,S,'.b',t,S1,'-m',t,S4,'-c');
legend(Y,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_S(t)'

figure(15)
T=plot(t,E,'.b',t,E1,'-m',t,E4,'-c');
legend(T,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_E (t)'

figure(16)
f=plot(t,I,'.b',t,I1,'-m',t,I4,'-c');
legend(f,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_I(t)'

figure(17)
C=plot(t,C,'.b',t,C1,'-m',t,C4,'-c');
legend(C,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_Q(t)'



figure(18)
h=plot(t,R,'.b',t,R1,'-m',t,R4,'-c');
legend(h,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'H_R(t)'

figure(19)
H=plot(t,D,'.b',t,D1,'-m',t,D4,'-k');
legend(H,'u_1=0,u_2=0','u_1 \neq 0,u_2=0','u_2 \neq 0,u_1=0')
xlabel 'Time (days)'
ylabel 'D(t)'


% calculates the differential rates used in the integration.
function dpop=Diff_2_6(t,pop, parameter)
neta=parameter(1);omega=parameter(2);gamma=parameter(3);beta=parameter(4);mu=parameter(5);
delta=parameter(6);rho=parameter(7);alpha=parameter(8);N=parameter(9);u=parameter(10);
u3=parameter(11);
S=pop(1);E=pop(2);I=pop(3);C=pop(4);D=pop(5);R=pop(6);S1=pop(7);E1=pop(8);I1=pop(9);
C1=pop(10);D1=pop(11);R1=pop(12);S4=pop(13);E4=pop(14);I4=pop(15);C4=pop(16);D4=pop(17);R4=pop(18);
dpop=zeros(18,1);
dpop(1)=omega-(beta*I*S)./(N)-mu*S;
dpop(2)=(beta*(I)*S)./(N)-(neta+mu)*E;
dpop(3)=neta*E-(gamma+mu+delta+alpha).*I;
dpop(4)=alpha.*I-(rho+mu+delta).*C;
dpop(5)=delta*(I+C);
dpop(6)=rho*C+gamma*I-mu*R;
dpop(7)=omega-(1-u).*((beta*I1*S1)./(N))-mu*S1;
dpop(8)=(1-u).*((beta*I1*S1)./(N))-(neta+mu)*E1;
dpop(9)=neta*E1-(gamma+mu+delta+alpha).*I1;
dpop(10)=alpha*I1-(rho+mu+delta).*C1;
dpop(11)=delta*(I1+C1);
dpop(12)=rho*C1+gamma*I1-mu*R1;


dpop(13)=omega-((beta*I4*S4)./(N))-(mu+u3)*S4;
dpop(14)=((beta*I4*S4)./(N))-(neta+mu)*E4;
dpop(15)=neta*E4-(gamma+mu+delta+alpha).*I4;
dpop(16)=alpha*I4-(rho+mu+delta).*C4;
dpop(17)=delta*(I4+C4);
dpop(18)=rho*C4+gamma*I4-mu*R4+u3*S4;

