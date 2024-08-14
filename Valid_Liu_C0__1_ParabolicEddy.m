 
                                  %%  Code for Fig. 8(b) for parabolic eddy %% 

                    %%     Validation with Liu (2016) with Uniform sediment at the inlet  %%  
                               
clc; clear all; 

m = 18; %m+1 = Total of number of considered Chebyshev polynomial

                           %% Discretized data of Liu (2016) %% 

CC1P = [0.002303925,0.0741697,0.1625158,0.23636743,0.28912392,0.37618554,0.43291286,0.5001986,0.5529803,0.5978551,0.6532788,...
    0.68627423,0.7285175,0.78130877,0.83674896,0.8882227,0.9383702,0.989837,1.0610831,1.1455114,1.2272907,1.3156508,1.4593859,1.6123365];

ZZ1P = [0.9991189,0.9921042,0.9780462,0.95605224,0.9349301,0.88651025,0.84953076,0.8011022,0.7544294,0.7042288,0.65227085,...
    0.616162,0.56067395,0.5043095,0.43561134,0.37836525,0.32904807,0.27885044,0.22513726,0.17319202,0.13358045,0.10542551,0.074651785,0.057098005];
plot(CC1P,ZZ1P,'Marker','square','LineStyle','none','markersize',5,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta',...
    'DisplayName', 'Data of Liu (2016) for parabolic eddy');
hold on


                         %%     Generalized Value   %%    
Z0 = 0.001; %Non-dimensional start elevation
N=2000; %Number of grid points along stream-wise direction
dx=0.001; %Step size  
N*dx;
V0=0.2; %Non-dimensional Clear water settling velocity 
B=0.2; %Non-dimensional deposition velocity 
Cstar=1; %Non-dimensional equilibrium sediment concentration
kappa=0.41; %von-karman constant
za=0.05; %Non-dimensional reference height
C0 = 1; %Inlet concentration
n = 50; %number of grid points along the vertical direction
dz=1/n; % Step size in the vertical  

US = zeros(m+1, m+1);
US1= zeros(m+1,m+1);
DUS1=zeros(m+1,m+1);
D2US1=zeros(m+1,m+1);

US(1,1)=1;
for i = 2:m+1
    for k1 =0:i-1
        US(i,k1+1) = (-1)^k1 * 2^(2*(i-1)-2*k1) * factorial(2*(i-1)-k1+1)/(factorial(k1) * factorial(2*(i-1)-2*k1+1));
    end
end

XR = roots(US(m,1:m));
% XR'

US1(1,1) = 1;
for i = 2:m+1
    for k1 =0:i-1
        US1(i,k1+1) = (-1)^(i-1-k1) * 2^(2*k1) * factorial((i-1)+k1+1)/(factorial((i-1)-k1) * factorial(2*k1+1));
    end
end
% US1

for i = 1:m+1
    for k1 =1:i-1
        DUS1(i,k1+1) = (-1)^(i-1-k1) * 2^(2*k1) * k1 * factorial((i-1)+k1+1)/(factorial((i-1)-k1) * factorial(2*k1+1));
    end
end
% DUS1

for i = 1:m+1
    for k1 =2:i-1
        D2US1(i,k1+1) = (-1)^(i-1-k1) * 2^(2*k1) * k1*(k1-1) * factorial((i-1)+k1+1)/(factorial((i-1)-k1) * factorial(2*k1+1));
    end
end
% D2US1

K = zeros(m-1,m+1); 

for r = 1:m-1
    for j = 1:m+1
        if j == m
            K(r,j) = 0; 
        else
            for i = 1:m+1
                K(r,j) = K(r,j) + US1(j,i)*XR(r)^(i-1);
            end
        end
    end
end

R1 = zeros(m-1,m+1);
R2 = zeros(m-1,m+1);
R = zeros(m-1,m+1);
for r = 1:m-1
    for j = 1:m+1
        for i = 2:m+1
            R1(r,j) = R1(r,j) + DUS1(j,i)*XR(r)^(i-2); 
        end
    end
end

for r = 1:m-1
    for j = 1:m+1
        for i = 3:m+1
            R2(r,j) = R2(r,j) + D2US1(j,i)*XR(r)^(i-3);  
        end
    end
end


for r = 1:m-1        
    R(r,:) = (Z0-log(Z0)-1)/((1-Z0)*log((za+(1-za)*XR(r))/Z0))* kappa*(1-XR(r))*(za+(1-za)*XR(r))/(1-za)*R2(r,:) + ...
        (Z0-log(Z0)-1)/((1-Z0)*log((za+(1-za)*XR(r))/Z0))* ( kappa*(-za+(1-za)*(1-2*XR(r)) )/(1-za) + V0/(1-za) )*R1(r,:); 
end
PP = zeros(m+1,m+1);

for j = 1:m+1
    for i=2:m+1
        PP(m,j) = PP(m,j) + kappa*(1-1)*(za+(1-za)*1)*DUS1(j,i);
    end
    PP(m,j) = PP(m,j) + sum(US1(j,:))*V0;
end
    
for j=1:m+1
    PP(m+1,j) = PP(m+1,j) + kappa*(1-0)*(za+(1-za)*0)*DUS1(j,2) + (V0-B)*US1(j,1); 
end


                        %%        For matrix loop   %%

QQ = zeros(m+1,m+1);
bb=zeros(m+1,1);
CC=zeros(m+1,N+1);
CC(1,1) = C0; 

for nn=1:N
    TT(nn)=nn*dx;
    
    for ii = 1:m-1
        for jj = 1:m+1
            PP(ii,jj) = K(ii,jj) - dx*R(ii,jj);
            QQ(ii,jj) = K(ii,jj);
        end 
    end
    
    bb(m+1,:) = -B*Cstar; 
    
    P1Q=inv(PP)*QQ;
    CC(:,nn+1) = P1Q*CC(:,nn) + PP\bb;
end 

CC(:,N+1)'

SpecP1Q=max(abs(eig(P1Q)))

ZZ = zeros(n+1,1);
CCm = zeros(m+1,n+1);
CCv = zeros(n+1,1);

for r=1:n+1
    ZZ(r) = (r-1)*dz;
    for i=1:m+1
        for j = 1:m+1
            CCm(i,r) = CCm(i,r) + US1(i,j)*ZZ(r).^(j-1); 
        end
        CCv(r) = CCv(r) + CC(i,N+1).*CCm(i,r); 
    end
    ZZ(r) = za+(1-za)*ZZ(r);
end

ZZ'
% CCv'

plot(CCv,ZZ,'b','DisplayName', 'Present solution') 
hold on 


