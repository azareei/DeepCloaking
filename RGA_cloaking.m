%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%REAL-CODED GENETIC ALGORITHM
%%%%      combination with 
%%%%            UNDX (UNIMODAL NORMAL DISTRIBUTION CROSSOVER)
%%%%      and   MGG  (MINIMAL GENERATION GAP) 
%%%%                   coded by takahito Iida on 2017/12/19
%%%%                            final revision is 2021/6/10
%%%%                                          at Osaka University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% data set%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Information of calculation%%%%%%%%%%%%%%%
NI = 40;  % number of individuals, recommend: 15n:50n
NC = 12;   % number of cross over 2*NC+2, reccomend: 10n
NG = 6000; % number of generation
NL=4; %number of layer
NP = 2*NL;   % number of variables

%%%%%%%%%control parameters%%%%%%%%%%%%%%%
min_solution=1000000; %initial value of solution
min_parameter=zeros(1,NP); %initial value of parameters


low_limit_parameter=ones(1,NP); %lower limit value of parameters
range_parameter=ones(1,NP); %range of values of parameters

low_limit_parameter(1:NL)=0.01;  %gamma
range_parameter(1:NL)=0.5;          %gamma

low_limit_parameter(NL+1:2*NL)=0.01;  %beta
range_parameter(NL+1:2*NL)=0.5;          %beta


%%%%%%%%%gene initial values%%%%%%%%%%%%%%%
pgene=zeros(NI,NP);   %population
fgene=zeros(2*NC+2,NP); %parents+children


%% RGA calculation%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%initial population build
[pgene]=initial_population(NP,NI,low_limit_parameter,range_parameter,pgene);

%%%RGA start
for ii=1:NG
    number_generation=ii
    [selected]=copy_select(NI);
    [fgene]=family_product(NP,NC,selected,pgene,fgene);
    [pgene,min_solution,min_parameter]=evaluation(NP,NC,selected,pgene,fgene,min_solution,min_parameter);

    plot(ii,min_solution,'.'); hold on; 
end

'optimised parameters'
min_parameter
min_solution

%% initial population builds%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pgene]=initial_population(NP,NI,low_limit_parameter,range_parameter,pgene)


for jj=1:NI
    for ii=1:NP
       pgene(jj,ii)=rand*range_parameter(1,ii)+low_limit_parameter(1,ii);
    end
end



end %end for function
%% selection of copy individuals%%%%%%
%%%% three individuals are selected without replacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [selected]=copy_select(NI)
selected=zeros(1,3);
y=randperm(NI,3);

for jj=1:3
    selected(jj)=y(jj);
end

end %end for function
%% children generate%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fgene]=family_product(NP,NC,selected,pgene,fgene)

%%%parents selection
fgene(1,:)=pgene(selected(1),:); %x1
fgene(2,:)=pgene(selected(2),:); %x2
sgene(1,:)=pgene(selected(3),:); %x3

middle=zeros(1,NP);
dvector=zeros(1,NP);
x31=zeros(1,NP);
dx31=zeros(1,NP);
sub1=0;
sub2=0;
sub3=0;

%% UNIMODAL NORMAL DISTRIBUTION CROSSOVER
for ii=1:NP
    middle(ii)=(fgene(1,ii)+fgene(2,ii))/2;
    dvector(ii)=fgene(2,ii)-fgene(1,ii);
    x31(ii)=sgene(1,ii)-fgene(1,ii);
    dx31(ii)=dvector(ii)*x31(ii);
    sub1=sub1+x31(ii)^2;
    sub2=sub2+dx31(ii)^2;    
    sub3=sub3+dvector(ii)^2;    
end

dlength=sqrt(sub1-sub2/sub3); %D

%%%preparation for ORTHONORMAL BASIS
%%%save max value of d
maxd=0;
for ii=1:NP
    if maxd < abs(dvector(ii))
        maxd=abs(dvector(ii));
        maxnumber=ii;
    end
end
%%%preparation of NP independent vector
A=eye(NP); %identity matrix
A(maxnumber,:)=dvector;
AA=zeros(NP);
AA(1,:)=A(maxnumber,:);  %AA(1,:) is dvector
for ii=2:NP
    AA(ii,:)=(ii <= maxnumber)*A(ii-1,:)+(ii > maxnumber)*A(ii,:);
end


%%%GRAM SCHMIDT ORTHONORMALIZATION
f=AA; %f(k)=a(k)
for k=1:NP %e(k) ORTHONORMAL BASIS
    if k > 1 
        for ii=1:k-1
            f(k,:)=f(k,:)-dot(AA(k,:),e(ii,:))*e(ii,:); %f(k)=a(k)-sum(a(k)*e(ii))*e(ii)
        end % end for ii
    end % end for if k>1
    
    flength=0;
    for jj=1:NP % |f|
        flength=flength+f(k,jj)^2;
    end %end for jj
    
e(k,:)=f(k,:)/sqrt(flength); % ORTHONORMAL BASIS
end % end for k


%%%children production
sigma1=0.5; %alpha
sigma2=0.35/sqrt(NP); %beta

for ii=1:NC
    w1=sigma1*randn;
    w2=sigma2*randn;
    fgene(2+ii,:)=middle+w1*dvector;
    fgene(2+NC+ii,:)=middle-w1*dvector;
    for k=2:NP
        fgene(2+ii,:)=fgene(2+ii,:)+dlength*w2*e(k,:);
        fgene(2+NC+ii,:)=fgene(2+NC+ii,:)-dlength*w2*e(k,:);
    end
end

end %end for function

%% evaluation%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pgene,min_solution,min_parameter]=evaluation(NP,NC,selected,pgene,fgene,min_solution,min_parameter);

parameter=zeros(1,NP);
result=zeros(1,2+2*NC);
result2=zeros(1,2+2*NC);
abs_result=zeros(1,2+2*NC);
p=zeros(1,2+2*NC);
pp=zeros(1,2+2*NC);

min_result=1000; %initial value
max_result=0; %initial value

for nn=1:2+2*NC
    for ii=1:NP
        parameter(ii)=fgene(nn,ii);
    end
    
    judge=0; %judge for lethal gene
    [judge]=lethalgene(NP,parameter,judge);
    if judge > 5
        result(nn)=25;
    else
        [result(nn)]=function_cloaking(NP,parameter);
    end
    result2(nn)=result(nn)^2;
    abs_result(nn)=abs(result(nn));
    
    %%%minimum individual
    if abs_result(nn)< min_result
        min_result=abs_result(nn);
        min_number=nn;
    end
    %%%maximum individual
    if abs_result(nn)> max_result
        max_result=abs_result(nn);
        max_result2=result2(nn);
    end    
    %%%save minimum individual for total generation
    if abs_result(nn) < abs(min_solution)
        min_solution=result(nn);
        min_parameter=parameter;
    end
        
end
%%%fitness evaluation
sumfit=0;
for nn=1:2+2*NC
    sumfit=sumfit+(max_result2-result2(nn)); %summention of result
end
psave=0;
for nn=1:2+2*NC
    p(nn)=(max_result2-result2(nn))/sumfit; % probability
    psave=psave+p(nn);
    pp(nn)=psave; % stack probability
end

%%%ROULETTE SELECT
roulette=rand;
elite_number1=min_number;
elite_number2=~(elite_number1==1)*1+(elite_number1==1)*2; %initial value
for nn=1:2+2*NC
    if roulette < pp(nn) && nn ~= min_number
        elite_number2=nn;
    end
end

pgene(selected(1),:)=fgene(elite_number1,:);
pgene(selected(2),:)=fgene(elite_number2,:);

end %end for function
%% lethal gene%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [judge]=lethalgene(NP,parameter,judge)

for ii=1:NP/2 %gamma
    if parameter(ii)<=0.01
        judge=10;
    elseif parameter(ii)>0.5
        judge=10;
    end
end

for ii=NP/2+1:2*NP/2 %beta
    if parameter(ii)<=0.01
        judge=10;
    end
end


end %end for function
