clear all
close all
%% ABM
% In the agent based model we consider an initial amount of infected and non-infected macrophages, neutrophils and extracellular bacillus. 
euler_yes=false;
% We define the states of the different plots functions to know which graphs will be shown.
Plots_Nsimulations=false; Plots_abm=false; Plots_continuous = false; Plots_comparison = true; alpha_comp=false; plots3dcompa=false;

% We make N simulations of our ABM model to acquire enough data.
N = 500;
% We will save the data for later usage
BiMat=zeros(1000,N); MiMat=zeros(1000,N); LiMat=zeros(1000,N); 
MuMat=zeros(1000,N); NMat=zeros(1000,N); BeMat=zeros(1000,N);

% We define our initial conditions and execute the ibm code to realize
% the ABM model simulation. 
m_i0 =4; n0=7; m_u0=6; b_e0 =7; Tmax=100; dt=0.1; 
global delta
global gamm
global sigma
global eta
delta=30; gamm=0.01; sigma=20; eta=80; alph=60;
tic
[lisis,Bi,Mi,Mu,Be,Nvector,Ti] = ibm(m_i0,n0,m_u0,b_e0,Tmax,dt,delta,gamm,sigma,eta,alph);
toc
% We show the plots of the ABM model.
plotsabm(Plots_abm, Ti, Bi, Be, Mi, Mu, Nvector);

% We execute the N simulations and compute average results. 
for j=1:N
    % By calling the function ibm we make each simulation
   [lisis,Bi,Mi,Mu,Be,Nvector,Ti] = ibm(m_i0,n0,m_u0,b_e0,Tmax,dt,delta,gamm,sigma,eta,alph);
    BiMat(:,j)=Bi; MiMat(:,j)=Mi; LiMat(:,j)=(1/dt).*lisis; MuMat(:,j)=Mu;
    BeMat(:,j)=Be; NMat(:,j)=Nvector;
end

% Taking the mean value of our simulations:
BiA=sum(BiMat,2)./N; MiA=sum(MiMat,2)./N; LiA=sum(LiMat,2)./N; MuA=sum(MuMat,2)./N; BeA=sum(BeMat,2)./N; NA=sum(NMat,2)./N;
% We execute the plot function of the averaged results of N different simulations.
plotsN(Plots_Nsimulations, BiMat, MiMat, LiMat, BiA, MiA, LiA, N)

%% Continuous Model:
% In the continuous model we use the acquired data to adjust a continuous
% function to the previously discrete lysis term. And thus we relate
% effectively the number of lysis per day with bi and mi through the ABM
% model.
% We adjust the continuous function: f(x,y) = p00 + p10*x + p01*y + p20*x.^2 
% + p11*x.*y + p02*y.^2 + p21*x.^2.*y + p12*x.*y.^2 to the data.
[fitresult, ~] = fitting(BiA, MiA, LiA); 

% We define as global variables the coefficients of our fitting function, so we can
% use them in any function.
global p00
p00=fitresult.p00;
global p01
p01=fitresult.p01;
global p10
p10=fitresult.p10;
global p20
p20=fitresult.p20; 
global p11
p11=fitresult.p11; 
global p02
p02=fitresult.p02; 

% Initially the total number of lysis (L) is zero:
L=0; b_i0 = sum(randi(5,m_i0,1)); Tmax=100; dt = 0.1; yo=[b_i0, b_e0, m_u0, m_i0, n0, L]; 

% We integrate the ODE to find the evolution of our system:
if euler_yes
    [Bie, Bee, Mie, Mue, Nvectore, Tie]=euler([yo, Tmax, dt]);
else
    tic
    [Tie,Y] = ode15s(@odefun,[0 Tmax],yo);
    Bie=Y(:,1); Bee=Y(:,2); Mue=Y(:,3); Mie=Y(:,4); Nvectore=Y(:,5);
    toc
end

plotscontinuous(Plots_continuous, Tie, Bie, Bee, Mie, Mue, Nvectore);
plotscomparison(Plots_comparison,Ti,Tie,Bie,Bi,Bee,Be,Mue,Mu,Mie,Mi,Nvectore,NA);

%% Different alpha
% We repeat the same simulations changing the constant alpha in order to observe
% the effect on the results of our model (ABM (Diff_alph=true) and
% Continuous (Diff_alph=false)).
Diff_alph=false;
alph=30; N=500; dt=0.1; Tmax=100;
BiMat1=zeros(Tmax/dt,N); MiMat1=zeros(Tmax/dt,N); LiMat1=zeros(Tmax/dt,N); 
MuMat1=zeros(Tmax/dt,N); NMat1=zeros(Tmax/dt,N); BeMat1=zeros(Tmax/dt,N);

if Diff_alph == true
    % DISCRETE MODEL:
    % We start with the discrete model for alpha=30
    [lisis1,Bi1,Mi1,Mu1,Be1,Nvector1,Ti] = ibm(m_i0,n0,m_u0,b_e0,Tmax,dt,delta,gamm,sigma,eta,alph);
    Alphaplot(alpha_comp, Ti, Ti, Bi, Bi1, Be, Be1,Mi,Mi1, Mu, Mu1, Nvector, Nvector1 )
else
    %CONTINUOUS MODEL
    for j=1:N
        % By calling the function ibm we make each simulation
       [lisis1,Bi1,Mi1,Mu1,Be1,Nvector1,Ti] = ibm(m_i0,n0,m_u0,b_e0,Tmax,dt,delta,gamm,sigma,eta,alph);
        BiMat1(:,j)=Bi1; 
        MiMat1(:,j)=Mi1;
        LiMat1(:,j)=(1/dt).*lisis1; 
        MuMat1(:,j)=Mu1;
        BeMat1(:,j)=Be1; 
        NMat1(:,j)=Nvector1;
    end

    % Taking the mean value of our simulations:
    BiA1=sum(BiMat1,2)./N;
    MiA1=sum(MiMat1,2)./N; 
    LiA1=sum(LiMat1,2)./N;
    MuA1=sum(MuMat1,2)./N;
    BeA1=sum(BeMat1,2)./N;
    NA1=sum(NMat1,2)./N;

    %Then we calulate the continuous model
    [fitresult1, gof] = fitting(BiA1, MiA1, LiA1); 

    %We define again our new coefficients for the lysis fitting function.
    global p00a
    p00a=fitresult1.p00;
    global p01a
    p01a=fitresult1.p01;
    global p10a
    p10a=fitresult1.p10;
    global p20a
    p20a=fitresult1.p20; 
    global p11a
    p11a=fitresult1.p11; 
    global p02a
    p02a=fitresult1.p02; 

    % We execute the simulations of our continuous model for alpha=30
    L=0; b_i0 = sum(randi(5,m_i0,1));
    Tmax=100; yo1=[b_i0, b_e0, m_u0, m_i0, n0, L]; 

    [t1,Y1] = ode15s(@odefun_30,[0 Tmax],yo1);
    b_i1=Y1(:,1); b_e1=Y1(:,2); m_u1=Y1(:,3); m_i1=Y1(:,4); n1=Y1(:,5);
    Alphaplot(alpha_comp,Tie,t1,Bie,b_i1,Bee,b_e1,Mie,m_i1,Mue,m_u1,Nvectore,n1)
end
%% Comparation of the squared error R^2
% The fitting of the lysis function to the ABM data was made using only Bi
% and Mi as dependent variables. This choice was made after comparing the
% squared error of the adjustment of a plane to the IBM data considering
% the different possible dependences.
plotcomp3d(plots3dcompa,BiA,MiA,LiA,MuA,NA)

% Now we try to fitt the data that we have.
[~,gof_bm]=fitting_degree1(BiA, MiA, LiA);
[~,gof_mm]=fitting_degree1(MuA,MiA,LiA);
[~,gof_mn]=fitting_degree1(MuA,NA,LiA);
[~,gof_mb]=fitting_degree1(MuA,BiA,LiA);
[~,gof_nb]=fitting_degree1(NA,BiA,LiA);
[~,gof_nm]=fitting_degree1(NA,MiA,LiA);

R_bi_mi=gof_bm.rsquare; rmse_bi_mi=gof_bm.rmse; sse_bi_mi=gof_bm.sse; 
R_mu_mi=gof_mm.rsquare; rmse_mu_mi=gof_mm.rmse; sse_mu_mi=gof_mm.sse;
R_mu_n=gof_mn.rsquare; rmse_mu_n=gof_mn.rmse; sse_mu_n=gof_mn.sse;
R_mu_bi=gof_mb.rsquare; rmse_mu_bi=gof_mb.rmse; sse_mu_bi=gof_mb.sse;
R_n_bi=gof_nb.rsquare; rmse_n_bi=gof_nb.rmse; sse_n_bi=gof_nb.sse;
R_n_mi=gof_nm.rsquare; rmse_n_mi=gof_nm.rmse; sse_n_mi=gof_nm.sse;

RR=[R_bi_mi,R_mu_mi,R_mu_n,R_mu_bi,R_n_bi,R_n_mi]';
rmse=[rmse_bi_mi,rmse_mu_mi,rmse_mu_n,rmse_mu_bi,rmse_n_bi,rmse_n_mi]';
sse=[sse_bi_mi,sse_mu_mi,sse_mu_n,sse_mu_bi,sse_n_bi,sse_n_mi]';
names=["Bi_Mi","Mu_Mi","Mu_N","Mu_Bi","N_Bi","N_Mi"]';

%We compare two errors: R^2 and Root mean square error. The fitting is better 
% if R^2 is near 1 and RMSE & SSE as smaller the better.  
T=table(names,RR,rmse,sse)



%% Functions:

% ABM:
function [lisis,Bi,Mi,Mu,Be,Nvector,Ti] = ibm(m_i,n,m_u,b_e,Tmax,dt,delta,gamm,sigma,eta,alph)
    b_i = randi(5,m_i,1); u=log(2); ti=0.1;
    Bi=zeros(Tmax/dt,1); Mi=zeros(Tmax/dt,1); Mu=zeros(Tmax/dt,1); Be=zeros(Tmax/dt,1); Nvector=zeros(Tmax/dt,1);
    ii=1; Ti=zeros(Tmax/dt,1); 
    x=normrnd(alph,4,[length(b_i),1]); lisis=zeros(Tmax/dt,1);
    while ti<Tmax
        %Discrete EDOs:
        alpha=ceil(x);
        %Internal Bacillus Growth
        c_bi = u*b_i.*dt;%*(1 - (b_i./alpha).^3)*dt;
        %External Bacillus Growth
        c_be = u*b_e*(1 - (b_e./(delta*n)).^3)*dt; 
        %Uninfected Machrophages Phagocytosis
        fag_u = gamm*b_e*m_u*dt;
        %Infected Machrophages Phagocytosis
        fag_i = (b_i>0).*gamm.*b_e.*(1-(b_i./alpha).^3)*dt; 
        %Alveolar Saturation. It takes into account that the volume of an alveolus is finite. 
        g = 1 - ((length(b_i)+m_u)*4990 + n*299)/(5*10^6);
        %Machrophages Growth.
        flux_m = sigma*g*dt;
        %Neutrophiles Growth
        flux_n = eta*g*dt;

        %Conditions we must impose:
        %There cannot be more phagocytosed bacillus that the total ammount of bacillus present in the cell.
        fag_u = min(fag_u, m_u);
        fag_u = min(fag_u,b_e);

        % sum(fag_i) + fag_u < b_e 
        if (sum(fag_i) + fag_u) > b_e
            fag_i = (b_e - fag_u)/sum(fag_i>0)*(b_i>0);
        end
        
        %DISCRETIZATION OF OUR MODEL
        %Now we define a random vector of numbers between 1 and 0 that will
        %fix evcery step the probability of gaining 1 unit or remaining the same when
        %our number has decimals. 
        x1=rand(1,length(b_i));
        dbi = c_bi + fag_i;
        for i=1:length(b_i)
            if abs(dbi(i)-floor(dbi(i))) > x1(i)
            dbi(i)=floor(dbi(i))+1;
            else 
            dbi(i)=floor(dbi(i));
            end
        end
        %For the rest of variables, as they are just unique numbers, we
        %will impose the same conditions but now with a fixed number rand for every step.
        dbe = c_be - fag_u - sum(fag_i);
            if abs(dbe-floor(dbe)) > rand
            dbe=floor(dbe)+1;
            else 
            dbe=floor(dbe);
            end
        dmu = flux_m - fag_u;
            if abs(dmu-floor(dmu)) > rand
            dmu=floor(dmu)+1;
            else 
            dmu=floor(dmu);
            end
        dmi = fag_u;
            if abs(dmi-floor(dmi)) > rand
            dmi=floor(dmi)+1;
            else 
            dmi=floor(dmi);
            end
        dn = flux_n;
            if abs(dn-floor(dn)) > rand
            dn=floor(dn)+1;
            else 
            dn=floor(dn);
            end
        
        %We add the discrete variation:    
        m_i = m_i + dmi; Mi(ii)=m_i;
        b_i = b_i + dbi; Bi(ii)=sum(b_i);
        b_e = b_e + dbe; Be(ii)=b_e;
        m_u = m_u + dmu; Mu(ii)=m_u;
        n = n + dn; Nvector(ii)=n;

        %As dmi might not be an integer we wait until mi differs from the lenght of b_i an integer number using floor 
        RealLength = sum(b_i>0); 
        diff = m_i - RealLength; 

        if diff > 0
             %We create new infected macrophages:
              b_i=[b_i ; ones(floor(diff),1)];
        end

        %We kill the macrophagues that contain more than 60 bacillus
        x=[x ; normrnd(alph,4,[floor(diff),1])];
        %We create a new variables lysis that will save how many
        %machrophages are dying every step of time.
        lysis=0;
        for i=1:length(b_i)
            if b_i(i) > x(i)
                lysis=lysis+1;
                b_e = b_e + b_i(i);
                b_i(i)=0;
                m_i = m_i - 1; 
            end
        end
        %We update our lysis vector
        lisis(ii)=lysis; 
        %We update our time vector.
        Ti(ii) = ti;
        ti = ti + dt;
        ii=ii+1;
    end
end
function[]=plotsabm(Plots_abm, Ti, Bi, Be, Mi, Mu, Nvector)
%We execute the graphs of the abm model.
if Plots_abm
    figure(1)
    subplot(2,2,1)
    plot(Ti,Bi)
    grid on 
    xlabel('Time [days]')
    ylabel('Intracellular bacillus')

    subplot(2,2,2)
    plot(Ti,Be)
    grid on
    xlabel('Time [days]')
    ylabel('Extracellular bacillus')

    subplot(2,2,3)
    plot(Ti,Mu)
    grid on
    xlabel('Time [days]')
    ylabel('Un-infected macrophages')

    subplot(2,2,4)
    plot(Ti,Mi)
    grid on
    xlabel('Time [days]')
    ylabel('Infected macrophages')

    figure(2)
    plot(Ti,Nvector)
    grid on
    xlabel('Time [days]')
    ylabel('Neutrophils')
end
end
function[]=plotsN(Plots_Nsimulations, BiMat, MiMat, LiMat, BiA, MiA, LiA, N)
if Plots_Nsimulations
    %We plot the average of N different simulations of the ABM model, to
    %obtain smoother results.
    figure(3)
    plot3(BiMat,MiMat,LiMat,'o')
    grid on
    xlabel('Interior bacillus')
    ylabel('Infected macrophages')
    zlabel('Lysis')
    title([num2str(N),' Simulations'])
    
    figure(4)
    plot3(BiA,MiA,LiA,'o')
    grid on
    xlabel('Interior bacillus')
    ylabel('Infected macrophages')
    zlabel('Lysis')
    title(['Average of the ',num2str(N),' simulations'])
end
end

% Continuous model with alpha=60.
function[dy]=odefun(~,y)
    global delta
    global gamm
    global sigma
    global eta
    b_i=y(1); b_e=y(2); m_u=y(3); m_i=y(4); n=y(5); L = y(6); u=log(2); alpha = 60;
    %Intern Bacillus Growth
    c_bi = u*b_i.*(1 - (b_i./(m_i*alpha)).^3);
    %External Bacillus Growth
    c_be = u*b_e*(1 - (b_e./(delta*n)).^3); 
    %Uninfected Macrophages Phagocytosis
    fag_u = gamm*b_e*m_u;
    %Infectes Macrophages Phagocytosis
    fag_i = gamm.*b_e*m_i*(1-(b_i./(m_i*alpha)).^3); 
    %Alveolar Saturation (it takes into account that the volume of an alveolus is finite)
    g = 1 - ((m_i+L+m_u)*4990 + n*299)/(5*10^6);
    %Machrophages Growth
    flux_m = sigma*g;
    %Neutrophiles Growth
    flux_n = eta*g;
    %Now we compute our finite EDOs
    db_i = c_bi + fag_i + fag_u - max(alpha*lysfun(b_i,m_i),0);
    db_e = c_be - fag_u - fag_i + max(alpha*lysfun(b_i,m_i),0); 
    dm_u= flux_m - fag_u;
    dm_i = fag_u - max(lysfun(b_i,m_i),0);
    dn = flux_n;
    %We compute the number of lysis we have every step.
    dL = max(lysfun(b_i, m_i),0);
    %We define our new states after the time step has finished.
    dy=[db_i; db_e; dm_u; dm_i; dn ;dL ]; 
end
function[Bie, Bee, Mie, Mue, Nvectore, Tie]=euler(y)
    global delta
    global gamm
    global sigma
    global eta
    b_i=y(1); b_e=y(2); m_u=y(3); m_i=y(4); n=y(5); L=y(6); u=log(2); alpha = 60; dt=y(8); ii=1; Tmax=y(7);
    Bie=zeros(Tmax/dt,1); Mie=zeros(Tmax/dt,1); Mue=zeros(Tmax/dt,1); Bee=zeros(Tmax/dt,1); Nvectore=zeros(Tmax/dt,1); Tie=zeros(Tmax/dt,1);
    while ii*dt < Tmax+dt
        %Intern Bacillus Growth
        c_bi = u*b_i*(1 - (b_i/(m_i*alpha))^3);
        %External Bacillus Growth
        c_be = u*b_e*(1 - (b_e/(delta*n))^3); 
        %Uninfected Macrophages Phagocytosis
        fag_u = gamm*b_e*m_u;
        %Infected Macrophages Phagocytosis
        fag_i = gamm*b_e*m_i*(1-(b_i/(m_i*alpha))^3); 
        %Alveolar Saturation (it takes into account that the volume of an alveolus is finite)
        g = 1 - ((m_i+L+m_u)*4990 + n*299)/(5*10^6);
        %Machrophages Growth
        flux_m = sigma*g;
        %Neutrophiles Growth
        flux_n = eta*g;
        %Now we compute our finite EDOs
        db_i = (c_bi + fag_i + fag_u - max(alpha*lysfun(b_i,m_i),0))*dt;
                     if abs(db_i-floor(db_i)) > rand
                     db_i=floor(db_i)+1;
                     else 
                     db_i=floor(db_i);
                     end
        db_e = (c_be - fag_u - fag_i + max(alpha*lysfun(b_i,m_i),0))*dt; 
                     if abs(db_e-floor(db_e)) > rand
                     db_e=floor(db_e)+1;
                     else 
                     db_e=floor(db_e);
                     end
        dm_u= (flux_m - fag_u)*dt;
                     if abs(dm_u-floor(dm_u)) > rand
                     dm_u=floor(dm_u)+1;
                     else 
                     dm_u=floor(dm_u);
                     end
        dm_i = (fag_u - max(lysfun(b_i,m_i),0))*dt;
                     if abs(dm_i-floor(dm_i)) > rand
                     dm_i=floor(dm_i)+1;
                     else 
                     dm_i=floor(dm_i);
                     end
        dn = flux_n*dt;
                     if abs(dn-floor(dn)) > rand
                     dn=floor(dn)+1;
                     else 
                     dn=floor(dn);
                     end
        dL = max(lysfun(b_i,m_i),0)*dt;
        % We compute the number of lysis we have every step.
        % We add the discrete variation:    
        m_i = m_i + dm_i; Mie(ii)=m_i; 
        b_i = b_i + db_i; Bie(ii)=b_i;
        b_e = b_e + db_e; Bee(ii)=b_e;
        m_u = m_u + dm_u; Mue(ii)=m_u;
        n = n + dn; Nvectore(ii)=n;
        L = L + dL;

        Tie(ii)=ii*dt;
        ii=ii+1;
    end
end

% Lysis fitting function.
function [fitresult, gof] = fitting(BiA, MiA, LiA)
%CREATEFIT(BIA,MIA,LIA)
%We create out fitting function useing the preinstalled Matlab curve
%fitting application.
[xData, yData, zData] = prepareSurfaceData( BiA, MiA, LiA );

% Set up fittype and options.
ft = fittype( 'poly22' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );
end
function[lys]=lysfun(bi,mi)
x=bi; y=mi;
%We define again our global coefficients in order to use them in the
%fitting lysis function.
global p00 
global p01
global p10
global p20
global p11 
global p02
%We define our fitting fucntion.
lys=p00+p10*x+p01*y+p20*x.^2+p11*x.*y+p02*y.^2;
end
function[]=plotscontinuous(Plots_continuous, t, b_i, b_e, m_i, m_u, n)
if Plots_continuous
    figure(5)
    subplot(2,2,1)
    plot(t,b_i)
    grid on
    xlabel('Time [days]')
    ylabel('Intracellular bacillus')

    subplot(2,2,2)
    plot(t,b_e)
    grid on
    xlabel('Time [days]')
    ylabel('Extracellular bacillus')

    subplot(2,2,3)
    plot(t,m_u)
    grid on
    xlabel('Time [days]')
    ylabel('Un-infected macrophages')

    subplot(2,2,4)
    plot(t,m_i)
    grid on
    xlabel('Time [days]')
    ylabel('Infected macrophages')

    figure(6)
    plot(t,n)
    grid on
    xlabel('Time [days]')
    ylabel('Neutrophils')
end

end

% Continuous model with alpha=30.
function[dy]=odefun_30(~,y)
global delta
global gamm
global sigma
global eta
b_i=y(1); b_e=y(2); m_u=y(3); m_i=y(4); n=y(5); L = y(6); u=log(2); alpha = 30;
    %Intern Bacillus Growth
    c_bi = u*b_i.*(1 - (b_i./(m_i*alpha)).^3);
    %External Bacillus Growth
    c_be = u*b_e*(1 - (b_e./(delta*n)).^3); 
    %Uninfected Macrophages Phagocytosis
    fag_u = gamm*b_e*m_u;
    %Infectes Macrophages Phagocytosis
    fag_i = gamm.*b_e*m_i*(1-(b_i./(m_i*alpha)).^3); 
    %Alveolar Saturation (it takes into account that the volume of an alveolus is finite)
    g = 1 - ((m_i+L+m_u)*4990 + n*299)/(5*10^6);
    %Machrophages Growth
    flux_m = sigma*g;
    %Neutrophiles Growth
    flux_n = eta*g;
    %Now we compute our finite EDOs
db_i = c_bi + fag_i + fag_u - max(alpha*lysfun30(b_i,m_i),0);
db_e = c_be - fag_u - fag_i + max(alpha*lysfun30(b_i,m_i),0); 
dm_u= flux_m - fag_u;
dm_i = fag_u - max(lysfun30(b_i,m_i),0);
dn = flux_n;
%We compute the number of lysis we have every step.
dL = max(lysfun30(b_i, m_i),0);
%We define our new states after the time step has finished.
dy=[db_i; db_e; dm_u; dm_i; dn ;dL]; 
end

% We define a second Lysis fitting function for alpha=30
function[lys1]=lysfun30(bi,mi)
x=bi; y=mi; 
%We define again our global coefficients in order to use them in the
%fitting lysis function.
global p00a 
global p01a
global p10a
global p20a
global p11a 
global p02a

%We define our fitting fucntion.
lys1=p00a+p10a*x+p01a*y+p20a*x.^2+p11a*x.*y+p02a*y.^2;
end
function[]=Alphaplot(alpha_comp, t, t1, b_i, b_i1, b_e, b_e1, m_i, m_i1, m_u, m_u1, n, n1 )
if alpha_comp
    
    figure(9)
    subplot(2,2,1)
    plot(t,b_i)
    hold on
    plot(t1,b_i1)
    xlabel('Time [days]')
    ylabel('Intracellular bacillus')
    legend('\alpha=60','\alpha=30')
    grid on
    hold off
   

    subplot(2,2,2)
    plot(t,b_e)
    hold on
    plot(t1,b_e1)
    xlabel('Time [days]')
    ylabel('Extracellular bacillus')
    legend('\alpha=60','\alpha=30')
    grid on
    hold off

    subplot(2,2,3)
    plot(t,m_u)
    hold on
    plot(t1,m_u1)
    xlabel('Time [days]')
    ylabel('Un-infected macrophages')
    legend('\alpha=60','\alpha=30')
    grid on
    hold off

    subplot(2,2,4)
    plot(t,m_i)
    hold on
    plot(t1,m_i1)
    xlabel('Time [days]')
    ylabel('Infected macrophages')
    legend('\alpha=60','\alpha=30')
    grid on
    hold off

    figure(10)
    plot(t,n)
    hold on
    plot(t1,n1)
    xlabel('Time [days]')
    ylabel('Neutrophils')
    legend('\alpha=60','\alpha=30')
    grid on
    hold off

end
end

% Fitting with a plane
function [fitresult, gof] = fitting_degree1(BiA, MiA, LiA)
%CREATEFIT(BIA,MIA,LIA)
%We create out fitting function useing the preinstalled Matlab curve
%fitting application.
[xData, yData, zData] = prepareSurfaceData( BiA, MiA, LiA );

% Set up fittype and options.
ft = fittype( 'poly11' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );
end

% 3D Plots between all the variables
function []=plotcomp3d(plots3dcompa,BiA,MiA,LiA,MuA,NA)
if plots3dcompa
    figure (11)
    subplot(3,2,1)
    plot3(BiA,MiA,LiA,'o')
    xlabel('Intracellular bacillus')
    ylabel('Infected macrophages')
    zlabel('Lysis')
    grid on

    subplot(3,2,2)
    plot3(MuA,MiA,LiA,'o')
    xlabel('Un-infected macrophages')
    ylabel('Infected macrophages')
    zlabel('Lysis')
    grid on

    subplot(3,2,3)
    plot3(MuA,NA,LiA,'o')
    xlabel('Un-infected macrophages')
    ylabel('Neutrophils')
    zlabel('Lysis')
    grid on

    subplot(3,2,4)
    plot3(MuA,BiA,LiA,'o')
    xlabel('Un-infected macrophages')
    ylabel('Intracellular bacillus')
    zlabel('Lysis')
    grid on

    subplot(3,2,5)
    plot3(NA,BiA,LiA,'o')
    xlabel('Neutrophils')
    ylabel('Intracellular bacillus')
    zlabel('Lysis')
    grid on

    subplot(3,2,6)
    plot3(NA,MiA,LiA,'o')
    xlabel('Neutrophils')
    ylabel('Infected macrophages')
    zlabel('Lysis')
    grid on
end
end 

% Comparison between both models
function[]=plotscomparison(Plots_comparison, Ti,t,b_i,BiA,b_e,BeA,m_u, MuA, m_i, MiA, n, NA)
%We plot the comparison between the both models realized, just to know how
%close is out continuous model of our ABM model.
if Plots_comparison
    figure(7)
    subplot(2,2,1)
    plot(t,b_i)
    hold on
    plot(Ti,BiA)
    xlabel('Time [days]')
    ylabel('Intracellular bacillus')
    legend('Continuous','ABM')
    grid on
    hold off

    subplot(2,2,2)
    plot(t,b_e)
    hold on
    plot(Ti,BeA)
    xlabel('Time [days]')
    ylabel('Extracellular bacillus')
    legend('Continuous','ABM')
    grid on
    hold off

    subplot(2,2,3)
    plot(t,m_u)
    hold on
    plot(Ti,MuA)
    xlabel('Time [days]')
    ylabel('Un-infected macrophages')
    legend('Continuous','ABM')
    grid on
    hold off

    subplot(2,2,4)
    plot(t,m_i)
    hold on
    plot(Ti,MiA)
    xlabel('Time [days]')
    ylabel('Infected macrophages')
    legend('Continuous','ABM')
    grid on
    hold off

    figure(8)
    plot(t,n)
    hold on
    plot(Ti,NA)
    xlabel('Time [days]')
    ylabel('Neutrophils')
    legend('Continuous','ABM')
    grid on
    hold off
end
end