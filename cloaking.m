clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculation of wave field around a cylinder with a multi-layered carpet
%%% carpet is designed to cloak the cylinder at k0=1.0
%%% 
%%%                                coded by 
%%%                                       Takahito Iida, Osaka University
%%%                                       Ahmad Zareei,Harvard University
%%%                                 supervised by Reza Alam, UC. Berkeley
%%%                                 final version is coded on 2021 June 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% data set%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Information of calculation%%%%%%%%%%%%%%%
NN=15; % number of wave modes in r direction
NM=20; % number of modes in theta direction
NMI=NM; % number of modes in theta direction for incident wave

NI=4; % number of layers MAX=NT
NR=NI+1; % number of radius

%%%%%%%%%Information of Waves%%%%%%%%%%%%%%%
k_inc=1.0; %wave number k
lambda = 2*pi/k_inc; %wave length lambda=2pi/k

H = lambda; % water depth
alpha = k_inc*tanh(k_inc*H); % alpha=K=omega^2/g=ktanh(kh)

%%%%%%%%%Information of cylinder%%%%%%%%%%%%%%%
ra = 1; % innter radius a
rb = 5; % outer radius

%%%%%%%%%Information of elastic plate%%%%%%%%%%%%%%%
beta0=1.0; % fundamental flexual rigidity = beta 
nu=0.25*ones(1,NI); % Poisson ratio 
gamma = zeros(1,NI); % gamma is mass rho_s*t/(rho_w*L)
beta=ones(1,NI).*beta0;

% % %%%%%%N=4%%%%%
 gamma(1)=0.24973416144981400;
 gamma(2)=gamma(1);
 gamma(3)=gamma(1);
 gamma(4)=gamma(1);
  
 beta(1)=0.08847287945600070;
 beta(2)=0.01216240715643230;
 beta(3)=0.32347070267820000;
 beta(4)=0.01000000005416440;

% xlswrite('gamma.xlsx',gamma(:)) 
% xlswrite('beta.xlsx',beta(:)) 

%%%%%%%%%nondimensional values%%%%%%%%%%%%%%%
ka=k_inc*ra ; %nondimensional wave number ka
lambda_a=lambda/ra; %nondimensional wavelength

%% material properties for cloaking%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%layer number i starts from outer to innter

RR = linspace(rb,ra,NR); % R(i) outer radius of i layer
ssign = ((-1).^(1:NR)+1)/2; %ssign(i)=0 (odd number) or 1 (even number)


%% wave number%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%wave number for water waves%%%%%%%%%%%%%%%
[k0,kI] = disper_water_ND(H,NN,alpha);
k0 = double(k0);
k = [k0,kI]; %k=[k(0),k(1),k(2)...k(NN)]

%%%%%%%%%wave number for plate mu(n,i) for i-th layer%%%%%%%%%%%%%%%
mu = zeros(NN+1+2,NI);
for ii = 1:NI
    [tem0,temI,tem1,tem2,tem3,tem4] = disper_plate_with_mass_ND(H,NN, alpha, beta(ii), gamma(ii));
    tem0 = double(tem0); %+Real
    tem1 = double(tem1); %-Real+Imag
    tem2 = double(tem2); %-Real-Imag
    tem3 = double(tem3); %+Real+Imag
    tem4 = double(tem4); %+Real-Imag
    mu(:,ii) = [tem0;temI.';tem3;tem4]; %mu(n,i)=[mu(0,i),mu(1,i),mu(2,i)...mu(NN,i),mu(-1,i),mu(-2,i)]
    test=alpha/((1-alpha*gamma(ii)+beta(ii)*tem3^4))+tem3*tan(tem3*H);
end

% xlswrite('k.xlsx',k) 
mu_imag(:,:)=imag(mu(:,:));
% xlswrite('mu_real.xlsx',mu(:,:))
% xlswrite('mu_imag.xlsx',mu_imag(:,:))

%%%%%%%%C0%%%%%%%%%%%%%%%%%%%
C0=k0^2/(alpha+(k0^2-alpha^2)*H);

%% setting for matrix%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%coefficient for wave elevation 1/(beta(i)*mu(n,i)^4+1-alpha*gamma)
GG = @(muu,ii) 1./(beta(ii)*muu.^4 + 1 - alpha*gamma(ii));

%%%derivative of bessel functions
dH1 = @(n,x) 0.5*(besselh(n-1,1,x) - besselh(n+1,1,x));
dJ = @(n,x) 0.5*(besselj(n-1,x) - besselj(n+1,x));
dY = @(n,x) 0.5*(bessely(n-1,x) - bessely(n+1,x));
dI = @(n,x) 0.5*(besseli(n-1,x) + besseli(n+1,x));
dK = @(n,x) (-0.5)*(besselk(n-1,x) + besselk(n+1,x));

%%%x=[a,b,c]
anm=zeros(NN+1,2*NM+1);
bnm=zeros(NN+3,2*NM+1);
cnm=zeros(NN+3,2*NM+1);

%% Z-solutions%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Z-solution for water waves 
f = @(z,n) (n==1)*cosh(k(n)*(z+H))/cosh(k(n)*H)...
    +~(n==1)*cos(k(n)*(z+H));

%%%Z-solution for plate
F = @(z,n,i) (n==1)*cosh(mu(n,i)*(z+H))/cosh(mu(n,i)*H)...
    +(n== NN+2 | n== NN+3)*cos(mu(n,i)*(z+H))/cos(mu(n,i)*H)...
    +~(n==1 | n== NN+2 | n== NN+3)*cos(mu(n,i)*(z+H));

%%%integral of f and F
%%%Aln(l,n)=int(-H,0)f(z,l)*f(z,n)
Aln = zeros(NN+1,NN+1);
for ll = 1:NN+1
    for nn=1:NN+1
        kn=1i*k(nn)*(nn==1)+k(nn)*~(nn==1);
        Aln(ll,nn)=(nn==ll)*(sin(kn*H)*cos(kn*H)+kn*H)/(2*kn)...
            +~(nn==ll)*0;
        diva=(nn==1)*cos(kn*H)^2+~(nn==1)*1;
        Aln(ll,nn)=Aln(ll,nn)/diva;
    end
end

%%%Blni(l,n,i)=int(-H,0)f(z,l)*F(z,n,i)
Blni = zeros(NN+1,NN+3,NI);
for ii = 1:NI
    for ll = 1:NN+1
    kl=(ll==1)*1i*k(ll)+~(ll==1)*k(ll);
    divb1=(ll==1)*cos(kl*H)+~(ll==1)*1;
        for nn=1:NN+3
            mun=1i*mu(nn,ii)*(nn==1)+mu(nn,ii)*~(nn==1);
            Blni(ll,nn,ii)=(kl*sin(kl*H)*cos(mun*H)-mun*cos(kl*H)*sin(mun*H))/(kl^2-mun^2);
            divb2=(nn==1 | nn== NN+2 | nn== NN+3)*cos(mun*H)+~(nn==1 | nn== NN+2 | nn== NN+3)*1;
            Blni(ll,nn,ii)=Blni(ll,nn,ii)/(divb1*divb2);
        end
    end
end

%% matrix building and solve it%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_unknown = NN+1+2*NI*(NN+3) %N+1+2*NT*(N+3)
coef = zeros(num_unknown,2*NM+1);  % x in Ax=b

for mm=-NM:NM
%for mm=1:1    
    'now calculating the m number'  
    mm
    %Ax=b
    A = zeros(num_unknown,num_unknown);
    b=zeros(num_unknown,1);

%%%%%%Continuity of phi - outer layer%%%%%%%%%
    cao=zeros(1,NN+1); %coefficient of a(n,m)
    cbo=zeros(1,NN+3); %coefficient of b(n,m)
    cco=zeros(1,NN+3); %coefficient of c(n,m)

    for nn =1:NN+1 
        cao(nn) = (nn==1)*besselh(mm,1,k(nn)*RR(1))/besselh(mm,1,k(nn)*RR(1))...
            +~(nn==1)*besselk(mm,k(nn)*RR(1))/besselk(mm,k(nn)*RR(1));                 
    end

    for nn = 1:NN+3
        cbo(nn) = (nn==1)*besselj(mm,mu(nn,1)*RR(1))/besselj(mm,mu(nn,1)*RR(1))...
            +~(nn==1)*besseli(mm,mu(nn,1)*RR(1))/besseli(mm,mu(nn,1)*RR(1));       
        
        cco(nn) = (nn==1)*besselh(mm,1,mu(nn,1)*RR(1))/besselh(mm,1,mu(nn,1)*RR(2))...
            +~(nn==1)*besselk(mm,mu(nn,1)*RR(1))/besselk(mm,mu(nn,1)*RR(2));
    end    

    
    for ll=1:NN+1
        for jj=1:NN+1+NN+3+NN+3
            if jj <= NN+1       
                A(ll,jj)=Aln(ll,jj)*cao(jj);
            elseif jj <= NN+1+NN+3
                j2=jj-(NN+1);
                A(ll,jj)=-Blni(ll,j2,1)*cbo(j2);
            else
                j2=jj-(NN+1+NN+3);
                A(ll,jj)=-Blni(ll,j2,1)*cco(j2);
            end
        end        
    end

    for ll=1:NN+1
        b(ll,1)=-(1i)^mm*besselj(mm,k(1)*rb)*Aln(ll,1);
    end
%%%%%%Continuity of dphi/dr - outer most layer%%%%%%%%%    
    cao_r=zeros(1,NN+1); %coefficient of a(n,m)
    cbo_r=zeros(1,NN+3); %coefficient of b(n,m)
    cco_r=zeros(1,NN+3); %coefficient of c(n,m)    
    
    for nn =1:NN+1
        cao_r(nn) = (nn==1)*k(nn)*dH1(mm,k(nn)*RR(1))/besselh(mm,1,k(nn)*RR(1))...
            +~(nn==1)*k(nn)*dK(mm,k(nn)*RR(1))/besselk(mm,k(nn)*RR(1));        
    end       

    for nn = 1:NN+3        
        cbo_r(nn) = (nn==1)*mu(nn,1)*dJ(mm,mu(nn,1)*RR(1))/besselj(mm,mu(nn,1)*RR(1))...
            +~(nn==1)*mu(nn,1)*dI(mm,mu(nn,1)*RR(1))/besseli(mm,mu(nn,1)*RR(1));
        
        cco_r(nn) = (nn==1)*mu(nn,1)*dH1(mm,mu(nn,1)*RR(1))/besselh(mm,1,mu(nn,1)*RR(2))...
            +~(nn==1)*mu(nn,1)*dK(mm,mu(nn,1)*RR(1))/besselk(mm,mu(nn,1)*RR(2));
    end       
   
    for ll=NN+2:NN+1+NN+1
        nl=ll-(NN+1);
        for jj=1:NN+1+NN+3+NN+3
            if jj <= NN+1       
                A(ll,jj)=Aln(nl,jj)*cao_r(jj);
            elseif jj <= NN+1+NN+3
                j2=jj-(NN+1);
                A(ll,jj)=-Blni(nl,j2,1)*cbo_r(j2);
            else
                j2=jj-(NN+1+NN+3);
                A(ll,jj)=-Blni(nl,j2,1)*cco_r(j2);
            end
        end        
    end    

    for ll=NN+2:NN+1+NN+1
        nl=ll-(NN+1);
        b(ll,1)=-(1i)^mm*k(1)*dJ(mm,k(1)*RR(1))*Aln(nl,1);
    end 
    
%%%%%%zero moment M at the outer most boundary%%%%%%%%%    
    cbom=zeros(1,NN+3); %coefficient of b(n,m)
    ccom=zeros(1,NN+3); %coefficient of c(n,m)

    for nn = 1:NN+3
        cbom(nn) = beta(1)/beta0*F(0,nn,1)*GG(mu(nn,1),1)*(...
            (nn==1)*(mu(nn,1)^2*besselj(mm,mu(nn,1)*RR(1))...
            -(1-nu(1))/RR(1)*(mu(nn,1)*dJ(mm,mu(nn,1)*RR(1))-mm^2/RR(1)*besselj(mm,mu(nn,1)*RR(1))))/besselj(mm,mu(nn,1)*RR(1))...
            +~(nn==1)*(mu(nn,1)^2*besseli(mm,mu(nn,1)*RR(1))...
            -(1-nu(1))/RR(1)*(mu(nn,1)*dI(mm,mu(nn,1)*RR(1))-mm^2/RR(1)*besseli(mm,mu(nn,1)*RR(1))))/besseli(mm,mu(nn,1)*RR(1)));        

        ccom(nn) = beta(1)/beta0*F(0,nn,1)*GG(mu(nn,1),1)*(...
            (nn==1)*(mu(nn,1)^2*besselh(mm,1,mu(nn,1)*RR(1))...
            -(1-nu(1))/RR(1)*(mu(nn,1)*dH1(mm,mu(nn,1)*RR(1))-mm^2/RR(1)*besselh(mm,1,mu(nn,1)*RR(1))))/besselh(mm,1,mu(nn,1)*RR(2))...
            +~(nn==1)*(mu(nn,1)^2*besselk(mm,mu(nn,1)*RR(1))...
            -(1-nu(1))/RR(1)*(mu(nn,1)*dK(mm,mu(nn,1)*RR(1))-mm^2/RR(1)*besselk(mm,mu(nn,1)*RR(1))))/besselk(mm,mu(nn,1)*RR(2)));
    
    end   
  
    
    lm=2*(NN+1)+1;
    for jj=NN+1+1:NN+1+2*(NN+3)
        j2=jj-(NN+1);
        if j2 <= NN+3
            A(lm,jj)=cbom(j2);
        else
            j3=j2-(NN+3);                
            A(lm,jj)=ccom(j3);
        end
    end              

%%%%%%zero shear force V at the outer most boundary%%%%%%%%%    
    cbov=zeros(1,NN+3); %coefficient of b(n,m)
    ccov=zeros(1,NN+3); %coefficient of c(n,m)    

    for nn = 1:NN+3
        cbov(nn) = beta(1)/beta0*F(0,nn,1)*GG(mu(nn,1),1)*(...
            (nn==1)*(mu(nn,1)^3*dJ(mm,mu(nn,1)*RR(1))...
            +(1-nu(1))/(RR(1)^2)*mm^2*(-mu(nn,1)*dJ(mm,mu(nn,1)*RR(1))+1/RR(1)*besselj(mm,mu(nn,1)*RR(1))))/besselj(mm,mu(nn,1)*RR(1))...
            +~(nn==1)*(mu(nn,1)^3*dI(mm,mu(nn,1)*RR(1))...
            +(1-nu(1))/(RR(1)^2)*mm^2*(-mu(nn,1)*dI(mm,mu(nn,1)*RR(1))+1/RR(1)*besseli(mm,mu(nn,1)*RR(1))))/besseli(mm,mu(nn,1)*RR(1)));
        
        ccov(nn) = beta(1)/beta0*F(0,nn,1)*GG(mu(nn,1),1)*(...
            (nn==1)*(mu(nn,1)^3*dH1(mm,mu(nn,1)*RR(1))...
            +(1-nu(1))/(RR(1)^2)*mm^2*(-mu(nn,1)*dH1(mm,mu(nn,1)*RR(1))+1/RR(1)*besselh(mm,1,mu(nn,1)*RR(1))))/besselh(mm,1,mu(nn,1)*RR(2))...
            +~(nn==1)*(mu(nn,1)^3*dK(mm,mu(nn,1)*RR(1))...
            +(1-nu(1))/(RR(1)^2)*mm^2*(-mu(nn,1)*dK(mm,mu(nn,1)*RR(1))+1/RR(1)*besselk(mm,mu(nn,1)*RR(1))))/besselk(mm,mu(nn,1)*RR(2)));   
    end
    
    
     lv=2*(NN+1)+1+1;
    for jj=NN+1+1:NN+1+2*(NN+3)
        j2=jj-(NN+1);            
        if j2 <= NN+3
            A(lv,jj)=cbov(j2);
        else
            j3=j2-(NN+3);                
            A(lv,jj)=ccov(j3);
        end
    end            

%%%%%%Continuity of phi - between layers%%%%%%%%%

    for ii=1:NI-1
    cbith=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    ccith=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbi1th=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    cci1th=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3
            cbith(nn) = (nn==1)*besselj(mm,mu(nn,ii)*RR(ii+1))/besselj(mm,mu(nn,ii)*RR(ii))...
            +~(nn==1)*besseli(mm,mu(nn,ii)*RR(ii+1))/besseli(mm,mu(nn,ii)*RR(ii));       
            
            ccith(nn) = (nn==1)*besselh(mm,1,mu(nn,ii)*RR(ii+1))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
            +~(nn==1)*besselk(mm,mu(nn,ii)*RR(ii+1))/besselk(mm,mu(nn,ii)*RR(ii+1));            

            cbi1th(nn) = (nn==1)*besselj(mm,mu(nn,ii+1)*RR(ii+1))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
            +~(nn==1)*besseli(mm,mu(nn,ii+1)*RR(ii+1))/besseli(mm,mu(nn,ii+1)*RR(ii+1)); 

            cci1th(nn) = (nn==1)*besselh(mm,1,mu(nn,ii+1)*RR(ii+1))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
            +~(nn==1)*besselk(mm,mu(nn,ii+1)*RR(ii+1))/besselk(mm,mu(nn,ii+1)*RR(ii+2));        
        end

        for ll=2*(NN+1)+2+(ii-1)*(NN+1)+1:2*(NN+1)+2+ii*(NN+1)
            nl=ll-(2*(NN+1)+2+(ii-1)*(NN+1));
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=Blni(nl,j2,ii)*cbith(j2);
                    A(ll,jj+2*(NN+3))=-Blni(nl,j2,ii+1)*cbi1th(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=Blni(nl,j3,ii)*ccith(j3);
                    A(ll,jj+2*(NN+3))=-Blni(nl,j3,ii+1)*cci1th(j3);                    
                end
            end
        end
    end % end for i-th    

 %%%%%%Continuity of dphi/dr - between layers%%%%%%%%%

    for ii=1:NI-1
    cbith_r=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    ccith_r=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbi1th_r=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    cci1th_r=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3
            cbith_r(nn) = (nn==1)*mu(nn,ii)*dJ(mm,mu(nn,ii)*RR(ii+1))/besselj(mm,mu(nn,ii)*RR(ii))...
            +~(nn==1)*mu(nn,ii)*dI(mm,mu(nn,ii)*RR(ii+1))/besseli(mm,mu(nn,ii)*RR(ii));

            ccith_r(nn) = (nn==1)*mu(nn,ii)*dH1(mm,mu(nn,ii)*RR(ii+1))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
            +~(nn==1)*mu(nn,ii)*dK(mm,mu(nn,ii)*RR(ii+1))/besselk(mm,mu(nn,ii)*RR(ii+1));            
                    
            cbi1th_r(nn) = (nn==1)*mu(nn,ii+1)*dJ(mm,mu(nn,ii+1)*RR(ii+1))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
            +~(nn==1)*mu(nn,ii+1)*dI(mm,mu(nn,ii+1)*RR(ii+1))/besseli(mm,mu(nn,ii+1)*RR(ii+1));
        
            cci1th_r(nn) = (nn==1)*mu(nn,ii+1)*dH1(mm,mu(nn,ii+1)*RR(ii+1))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
            +~(nn==1)*mu(nn,ii+1)*dK(mm,mu(nn,ii+1)*RR(ii+1))/besselk(mm,mu(nn,ii+1)*RR(ii+2));              
        end
        
        for ll=2*(NN+1)+2+(NI-1)*(NN+1)+(ii-1)*(NN+1)+1:2*(NN+1)+2+(NI-1)*(NN+1)+ii*(NN+1)
            nl=ll-(2*(NN+1)+2+(NI-1)*(NN+1)+(ii-1)*(NN+1));
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=Blni(nl,j2,ii)*cbith_r(j2);
                    A(ll,jj+2*(NN+3))=-Blni(nl,j2,ii+1)*cbi1th_r(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=Blni(nl,j3,ii)*ccith_r(j3);
                    A(ll,jj+2*(NN+3))=-Blni(nl,j3,ii+1)*cci1th_r(j3);                    
                end
            end
        end
    end % end for i-th   

    
%%%%%%Continuity of moment M - between layers%%%%%%%%%

    for ii=1:NI-1
    cbmith=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    ccmith=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbmi1th=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    ccmi1th=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3
            cbmith(nn) = beta(ii)/beta0*F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*(mu(nn,ii)^2*besselj(mm,mu(nn,ii)*RR(ii+1))...
                -(1-nu(ii))/RR(ii+1)*(mu(nn,ii)*dJ(mm,mu(nn,ii)*RR(ii+1))-mm^2/RR(ii+1)*besselj(mm,mu(nn,ii)*RR(ii+1))))/besselj(mm,mu(nn,ii)*RR(ii))...
                +~(nn==1)*(mu(nn,ii)^2*besseli(mm,mu(nn,ii)*RR(ii+1))...
                -(1-nu(ii))/RR(ii+1)*(mu(nn,ii)*dI(mm,mu(nn,ii)*RR(ii+1))-mm^2/RR(ii+1)*besseli(mm,mu(nn,ii)*RR(ii+1))))/besseli(mm,mu(nn,ii)*RR(ii)));        

            ccmith(nn) = beta(ii)/beta0*F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*(mu(nn,ii)^2*besselh(mm,1,mu(nn,ii)*RR(ii+1))...
                -(1-nu(ii))/RR(ii+1)*(mu(nn,ii)*dH1(mm,mu(nn,ii)*RR(ii+1))-mm^2/RR(ii+1)*besselh(mm,1,mu(nn,ii)*RR(ii+1))))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
                +~(nn==1)*(mu(nn,ii)^2*besselk(mm,mu(nn,ii)*RR(ii+1))...
                -(1-nu(ii))/RR(ii+1)*(mu(nn,ii)*dK(mm,mu(nn,ii)*RR(ii+1))-mm^2/RR(ii+1)*besselk(mm,mu(nn,ii)*RR(ii+1))))/besselk(mm,mu(nn,ii)*RR(ii+1)));
            
            cbmi1th(nn) = beta(ii+1)/beta0*F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*(mu(nn,ii+1)^2*besselj(mm,mu(nn,ii+1)*RR(ii+1))...
                -(1-nu(ii+1))/RR(ii+1)*(mu(nn,ii+1)*dJ(mm,mu(nn,ii+1)*RR(ii+1))-mm^2/RR(ii+1)*besselj(mm,mu(nn,ii+1)*RR(ii+1))))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
                +~(nn==1)*(mu(nn,ii+1)^2*besseli(mm,mu(nn,ii+1)*RR(ii+1))...
                -(1-nu(ii+1))/RR(ii+1)*(mu(nn,ii+1)*dI(mm,mu(nn,ii+1)*RR(ii+1))-mm^2/RR(ii+1)*besseli(mm,mu(nn,ii+1)*RR(ii+1))))/besseli(mm,mu(nn,ii+1)*RR(ii+1)));        

             ccmi1th(nn) = beta(ii+1)/beta0*F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*(mu(nn,ii+1)^2*besselh(mm,1,mu(nn,ii+1)*RR(ii+1))...
                -(1-nu(ii+1))/RR(ii+1)*(mu(nn,ii+1)*dH1(mm,mu(nn,ii+1)*RR(ii+1))-mm^2/RR(ii+1)*besselh(mm,1,mu(nn,ii+1)*RR(ii+1))))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
                +~(nn==1)*(mu(nn,ii+1)^2*besselk(mm,mu(nn,ii+1)*RR(ii+1))...
                -(1-nu(ii+1))/RR(ii+1)*(mu(nn,ii+1)*dK(mm,mu(nn,ii+1)*RR(ii+1))-mm^2/RR(ii+1)*besselk(mm,mu(nn,ii+1)*RR(ii+1))))/besselk(mm,mu(nn,ii+1)*RR(ii+2)));
        end                   

        ll=2*(NN+1)+2+2*(NI-1)*(NN+1)+ii;
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=cbmith(j2);
                    A(ll,jj+2*(NN+3))=-cbmi1th(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=ccmith(j3);
                    A(ll,jj+2*(NN+3))=-ccmi1th(j3);                    
                end
            end
    end % end for i-th

%%%%%%Continuity of shear force V - between layers%%%%%%%%%

    for ii=1:NI-1
    cbvith=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    ccvith=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbvi1th=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    ccvi1th=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3 
            cbvith(nn) = beta(ii)/beta0*F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*(mu(nn,ii)^3*dJ(mm,mu(nn,ii)*RR(ii+1))...
                +(1-nu(ii))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii)*dJ(mm,mu(nn,ii)*RR(ii+1))+1/RR(ii+1)*besselj(mm,mu(nn,ii)*RR(ii+1))))/besselj(mm,mu(nn,ii)*RR(ii))...
                +~(nn==1)*(mu(nn,ii)^3*dI(mm,mu(nn,ii)*RR(ii+1))...
                +(1-nu(ii))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii)*dI(mm,mu(nn,ii)*RR(ii+1))+1/RR(ii+1)*besseli(mm,mu(nn,ii)*RR(ii+1))))/besseli(mm,mu(nn,ii)*RR(ii)));
        
            ccvith(nn) = beta(ii)/beta0*F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*(mu(nn,ii)^3*dH1(mm,mu(nn,ii)*RR(ii+1))...
                +(1-nu(ii))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii)*dH1(mm,mu(nn,ii)*RR(ii+1))+1/RR(ii+1)*besselh(mm,1,mu(nn,ii)*RR(ii+1))))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
                +~(nn==1)*(mu(nn,ii)^3*dK(mm,mu(nn,ii)*RR(ii+1))...
                +(1-nu(ii))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii)*dK(mm,mu(nn,ii)*RR(ii+1))+1/RR(ii+1)*besselk(mm,mu(nn,ii)*RR(ii+1))))/besselk(mm,mu(nn,ii)*RR(ii+1)));
      
            cbvi1th(nn) = beta(ii+1)/beta0*F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*(mu(nn,ii+1)^3*dJ(mm,mu(nn,ii+1)*RR(ii+1))...
                +(1-nu(ii+1))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii+1)*dJ(mm,mu(nn,ii+1)*RR(ii+1))+1/RR(ii+1)*besselj(mm,mu(nn,ii+1)*RR(ii+1))))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
                +~(nn==1)*(mu(nn,ii+1)^3*dI(mm,mu(nn,ii+1)*RR(ii+1))...
                +(1-nu(ii+1))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii+1)*dI(mm,mu(nn,ii+1)*RR(ii+1))+1/RR(ii+1)*besseli(mm,mu(nn,ii+1)*RR(ii+1))))/besseli(mm,mu(nn,ii+1)*RR(ii+1)));
        
            ccvi1th(nn) = beta(ii+1)/beta0*F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*(mu(nn,ii+1)^3*dH1(mm,mu(nn,ii+1)*RR(ii+1))...
                +(1-nu(ii+1))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii+1)*dH1(mm,mu(nn,ii+1)*RR(ii+1))+1/RR(ii+1)*besselh(mm,1,mu(nn,ii+1)*RR(ii+1))))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
                +~(nn==1)*(mu(nn,ii+1)^3*dK(mm,mu(nn,ii+1)*RR(ii+1))...
                +(1-nu(ii+1))/(RR(ii+1)^2)*mm^2*(-mu(nn,ii+1)*dK(mm,mu(nn,ii+1)*RR(ii+1))+1/RR(ii+1)*besselk(mm,mu(nn,ii+1)*RR(ii+1))))/besselk(mm,mu(nn,ii+1)*RR(ii+2)));
        end

        ll=2*(NN+1)+2+2*(NI-1)*(NN+1)+(NI-1)+ii;
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=cbvith(j2);
                    A(ll,jj+2*(NN+3))=-cbvi1th(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=ccvith(j3);
                    A(ll,jj+2*(NN+3))=-ccvi1th(j3);                    
                end
            end
    end % end for i-th    
 
%%%%%%Continuity of wave elevation eta - between layers%%%%%%%%%

    for ii=1:NI-1
    cbeith=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    cceith=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbei1th=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    ccei1th=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3
            cbeith(nn) = F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*besselj(mm,mu(nn,ii)*RR(ii+1))/besselj(mm,mu(nn,ii)*RR(ii))...
                +~(nn==1)*besseli(mm,mu(nn,ii)*RR(ii+1))/besseli(mm,mu(nn,ii)*RR(ii)));  
              
            cceith(nn) = F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*besselh(mm,1,mu(nn,ii)*RR(ii+1))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
                +~(nn==1)*besselk(mm,mu(nn,ii)*RR(ii+1))/besselk(mm,mu(nn,ii)*RR(ii+1))); 
                                     
            cbei1th(nn) = F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*besselj(mm,mu(nn,ii+1)*RR(ii+1))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
                +~(nn==1)*besseli(mm,mu(nn,ii+1)*RR(ii+1))/besseli(mm,mu(nn,ii+1)*RR(ii+1))); 
                        
            ccei1th(nn) =  F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                 (nn==1)*besselh(mm,1,mu(nn,ii+1)*RR(ii+1))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
                 +~(nn==1)*besselk(mm,mu(nn,ii+1)*RR(ii+1))/besselk(mm,mu(nn,ii+1)*RR(ii+2)));                 
        end  

              
        ll=2*(NN+1)+2+2*(NI-1)*(NN+1)+2*(NI-1)+ii;
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=cbeith(j2);
                    A(ll,jj+2*(NN+3))=-cbei1th(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=cceith(j3);
                    A(ll,jj+2*(NN+3))=-ccei1th(j3);                    
                end
            end
    end % end for i-th        
    
%%%%%%Continuity of deta/dr - between layers%%%%%%%%%

    for ii=1:NI-1
    cbeith_r=zeros(1,NN+3); %coefficient of b(n,m,i) in ith layer
    cceith_r=zeros(1,NN+3); %coefficient of c(n,m,i)
    cbei1th_r=zeros(1,NN+3); %coefficient of b(n,m,i) in i+1th layer
    ccei1th_r=zeros(1,NN+3); %coefficient of c(n,m,i)
    
        for nn = 1:NN+3
            cbeith_r(nn) = F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                (nn==1)*mu(nn,ii)*dJ(mm,mu(nn,ii)*RR(ii+1))/besselj(mm,mu(nn,ii)*RR(ii))...
                +~(nn==1)*mu(nn,ii)*dI(mm,mu(nn,ii)*RR(ii+1))/besseli(mm,mu(nn,ii)*RR(ii)));
                        
            cceith_r(nn) = F(0,nn,ii)*GG(mu(nn,ii),ii)*(...
                 (nn==1)*mu(nn,ii)*dH1(mm,mu(nn,ii)*RR(ii+1))/besselh(mm,1,mu(nn,ii)*RR(ii+1))...
                 +~(nn==1)*mu(nn,ii)*dK(mm,mu(nn,ii)*RR(ii+1))/besselk(mm,mu(nn,ii)*RR(ii+1)));
        
            cbei1th_r(nn) = F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*mu(nn,ii+1)*dJ(mm,mu(nn,ii+1)*RR(ii+1))/besselj(mm,mu(nn,ii+1)*RR(ii+1))...
                +~(nn==1)*mu(nn,ii+1)*dI(mm,mu(nn,ii+1)*RR(ii+1))/besseli(mm,mu(nn,ii+1)*RR(ii+1)));
                         
            ccei1th_r(nn) =  F(0,nn,ii+1)*GG(mu(nn,ii+1),ii+1)*(...
                (nn==1)*mu(nn,ii+1)*dH1(mm,mu(nn,ii+1)*RR(ii+1))/besselh(mm,1,mu(nn,ii+1)*RR(ii+2))...
                +~(nn==1)*mu(nn,ii+1)*dK(mm,mu(nn,ii+1)*RR(ii+1))/besselk(mm,mu(nn,ii+1)*RR(ii+2))); 

        end  

        ll=2*(NN+1)+2+2*(NI-1)*(NN+1)+3*(NI-1)+ii;
            for jj=NN+1+(ii-1)*2*(NN+3)+1:NN+1+ii*2*(NN+3)
                j2=jj-(NN+1+(ii-1)*2*(NN+3));
                if j2 <= NN+3
                    A(ll,jj)=cbeith_r(j2);
                    A(ll,jj+2*(NN+3))=-cbei1th_r(j2);
                else
                    j3=j2-(NN+3);
                    A(ll,jj)=cceith_r(j3);
                    A(ll,jj+2*(NN+3))=-ccei1th_r(j3);                    
                end
            end
    end % end for i-th          
    
%%%%%%no flux condition at the inner most boundary dphi/dr=0%%%%%%%%%
    cbin_r=zeros(1,NN+3); %coefficient of b(n,m,i)
    ccin_r=zeros(1,NN+3); %coefficient of c(n,m,i)

    for nn = 1:NN+3
        cbin_r(nn) = (nn==1)*mu(nn,NI)*dJ(mm,mu(nn,NI)*RR(NI+1))/besselj(mm,mu(nn,NI)*RR(NI))...
            +~(nn==1)*mu(nn,NI)*dI(mm,mu(nn,NI)*RR(NI+1))/besseli(mm,mu(nn,NI)*RR(NI));
        
        ccin_r(nn) = (nn==1)*mu(nn,NI)*dH1(mm,mu(nn,NI)*RR(NI+1))/besselh(mm,1,mu(nn,NI)*RR(NI+1))...
            +~(nn==1)*mu(nn,NI)*dK(mm,mu(nn,NI)*RR(NI+1))/besselk(mm,mu(nn,NI)*RR(NI+1));     
    end        

    
    for ll=2*(NN+1)+2+2*(NN+1)*(NI-1)+4*(NI-1)+1:2*(NN+1)+2+2*(NN+1)*(NI-1)+4*(NI-1)+NN+1
        nl=ll-(2*(NN+1)+2+2*(NN+1)*(NI-1)+4*(NI-1));
        for jj=NN+1+(NI-1)*2*(NN+3)+1:NN+1+NI*2*(NN+3)
            j2=jj-(NN+1+(NI-1)*2*(NN+3));            
            if j2 <= NN+3
                A(ll,jj)=Blni(nl,j2,NI)*cbin_r(j2);
            else
                j3=j2-(NN+3);                
                A(ll,jj)=Blni(nl,j3,NI)*ccin_r(j3);
            end
        end        
    end

%%%%%%zero moment M at the inner most boundary%%%%%%%%%     
    cbim=zeros(1,NN+3); %coefficient of b(n,m,i)
    ccim=zeros(1,NN+3); %coefficient of c(n,m,i)
    
    for nn = 1:NN+3
        cbim(nn) = beta(NI)/beta0*F(0,nn,NI)*GG(mu(nn,NI),NI)*(...
            (nn==1)*(mu(nn,NI)^2*besselj(mm,mu(nn,NI)*RR(NI+1))...
            -(1-nu(NI))/RR(NI+1)*(mu(nn,NI)*dJ(mm,mu(nn,NI)*RR(NI+1))-mm^2/RR(NI+1)*besselj(mm,mu(nn,NI)*RR(NI+1))))/besselj(mm,mu(nn,NI)*RR(NI))...
            +~(nn==1)*(mu(nn,NI)^2*besseli(mm,mu(nn,NI)*RR(NI+1))...
            -(1-nu(NI))/RR(NI+1)*(mu(nn,NI)*dI(mm,mu(nn,NI)*RR(NI+1))-mm^2/RR(NI+1)*besseli(mm,mu(nn,NI)*RR(NI+1))))/besseli(mm,mu(nn,NI)*RR(NI)));
        
        ccim(nn) = beta(NI)/beta0*F(0,nn,NI)*GG(mu(nn,NI),NI)*(...
            (nn==1)*(mu(nn,NI)^2*besselh(mm,1,mu(nn,NI)*RR(NI+1))...
            -(1-nu(NI))/RR(NI+1)*(mu(nn,NI)*dH1(mm,mu(nn,NI)*RR(NI+1))-mm^2/RR(NI+1)*besselh(mm,1,mu(nn,NI)*RR(NI+1))))/besselh(mm,1,mu(nn,NI)*RR(NI+1))...
            +~(nn==1)*(mu(nn,NI)^2*besselk(mm,mu(nn,NI)*RR(NI+1))...
            -(1-nu(NI))/RR(NI+1)*(mu(nn,NI)*dK(mm,mu(nn,NI)*RR(NI+1))-mm^2/RR(NI+1)*besselk(mm,mu(nn,NI)*RR(NI+1))))/besselk(mm,mu(nn,NI)*RR(NI+1)));            
    end      

    
    lm=2*(NN+1)+2+2*(NN+1)*(NI-1)+4*(NI-1)+NN+1+1;
    for jj=NN+1+(NI-1)*2*(NN+3)+1:NN+1+NI*2*(NN+3)
        j2=jj-(NN+1+(NI-1)*2*(NN+3));            
        if j2 <= NN+3
            A(lm,jj)=cbim(j2);
        else
            j3=j2-(NN+3);                
            A(lm,jj)=ccim(j3);
        end
    end       
    
%%%%%%zero shear force V at the inner most boundary%%%%%%%%%     
    cbiv=zeros(1,NN+3); %coefficient of b(n,m,i)
    cciv=zeros(1,NN+3); %coefficient of c(n,m,i)

    for nn = 1:NN+3
        cbiv(nn) =  beta(NI)/beta0*F(0,nn,NI)*GG(mu(nn,NI),NI)*(...
            (nn==1)*(mu(nn,NI)^3*dJ(mm,mu(nn,NI)*RR(NI+1))...
            +(1-nu(NI))/(RR(NI+1)^2)*mm^2*(-mu(nn,NI)*dJ(mm,mu(nn,NI)*RR(NI+1))+1/RR(NI+1)*besselj(mm,mu(nn,NI)*RR(NI+1))))/besselj(mm,mu(nn,NI)*RR(NI))...
            +~(nn==1)*(mu(nn,NI)^3*dI(mm,mu(nn,NI)*RR(NI+1))...
            +(1-nu(NI))/(RR(NI+1)^2)*mm^2*(-mu(nn,NI)*dI(mm,mu(nn,NI)*RR(NI+1))+1/RR(NI+1)*besseli(mm,mu(nn,NI)*RR(NI+1))))/besseli(mm,mu(nn,NI)*RR(NI)));

        cciv(nn) = beta(NI)/beta0*F(0,nn,NI)*GG(mu(nn,NI),NI)*(...
            (nn==1)*(mu(nn,NI)^3*dH1(mm,mu(nn,NI)*RR(NI+1))...
            +(1-nu(NI))/(RR(NI+1)^2)*mm^2*(-mu(nn,NI)*dH1(mm,mu(nn,NI)*RR(NI+1))+1/RR(NI+1)*besselh(mm,1,mu(nn,NI)*RR(NI+1))))/besselh(mm,1,mu(nn,NI)*RR(NI+1))...
            +~(nn==1)*(mu(nn,NI)^3*dK(mm,mu(nn,NI)*RR(NI+1))...
            +(1-nu(NI))/(RR(NI+1)^2)*mm^2*(-mu(nn,NI)*dK(mm,mu(nn,NI)*RR(NI+1))+1/RR(NI+1)*besselk(mm,mu(nn,NI)*RR(NI+1))))/besselk(mm,mu(nn,NI)*RR(NI+1)));         
    end        

    
    lv=2*(NN+1)+2+2*(NN+1)*(NI-1)+4*(NI-1)+NN+1+1+1;
    for jj=NN+1+(NI-1)*2*(NN+3)+1:NN+1+NI*2*(NN+3)
        j2=jj-(NN+1+(NI-1)*2*(NN+3));            
        if j2 <= NN+3
            A(lv,jj)=cbiv(j2);
        else
            j3=j2-(NN+3);                
            A(lv,jj)=cciv(j3);
        end
    end    


    
%%%%%%solve matrix%%%%%%%%% 
%%%%%%Ax=b -> DAx=Db
    DD =  diag(1./max(abs(A.')));
    DA = DD*A;
    Db = DD*b;   
    
%%%%%%x=(DA)^(-1)*Db%%%%%%%
    coef(:,mm+NM+1) = DA\Db;    
        
end % end for m

for mm=1:2*NM+1
       for nn=1:NN+1
           anm(nn,mm)=coef(nn,mm);
       end       
       for ii=1:NI
           for nn=1:NN+3
               bnim(nn,ii,mm)=coef(nn+(ii-1)*2*(NN+3)+NN+1,mm);
               cnim(nn,ii,mm)=coef(nn+(ii-1)*2*(NN+3)+NN+3+NN+1,mm);
           end
       end              
end

%calculation of anm 
for mm=1:2*NM+1
       for nn=1:NN+1
           anm(nn,mm)=(nn==1)*coef(nn,mm)/besselh(mm-(NM+1),1,k(1)*rb)...
               +~(nn==1)*coef(nn,mm)/besselk(mm-(NM+1),k(nn)*rb);     
       end
       for ii=1:NI
           for nn=1:NN+3
               bnim(nn,ii,mm)=(nn==1)*coef(nn+(ii-1)*2*(NN+3)+NN+1,mm)/besselj(mm-(NM+1),mu(1,ii)*RR(ii))...
                   +~(nn==1)*coef(nn+(ii-1)*2*(NN+3)+NN+1,mm)/besseli(mm-(NM+1),mu(nn,ii)*RR(ii));
               cnim(nn,ii,mm)=(nn==1)*coef(nn+(ii-1)*2*(NN+3)+NN+3+NN+1,mm)/besselh(mm-(NM+1),1,mu(1,ii)*RR(ii+1))...
                   +~(nn==1)*coef(nn+(ii-1)*2*(NN+3)+NN+3+NN+1,mm)/besselk(mm-(NM+1),mu(nn,ii)*RR(ii+1));
           end
       end
end



%% scattered-wave energy at far field
S_energy=0;
T_energy=0;
I_energy=0;
for mm=1:2*NM+1
S_energy=S_energy+2*(anm(1,mm))*conj(anm(1,mm));
I_energy=I_energy+((1i)^(mm-(NM+1)))*conj(anm(1,mm))...
    +conj((1i)^(mm-(NM+1)))*(anm(1,mm));
end

S_energy=1/(2*C0*sqrt(alpha))*S_energy
I_energy=1/(2*C0*sqrt(alpha))*I_energy;
% I_energy
'Energy conservation'
T_energy=S_energy+I_energy



%% wave drift force at far field
F_drift=0;
for mm=-NM:NM-1   
    sub1=2*anm(1,mm+NM+1)*conj(anm(1,mm+NM+1+1));
    sub2=conj(1i^(mm+1))*anm(1,mm+NM+1);
    sub3=(1i^mm)*conj(anm(1,mm+NM+1+1));

    F_drift=F_drift+(sub1+sub2+sub3);
end
Fx_drift=k0/(C0*alpha)*imag(F_drift)



