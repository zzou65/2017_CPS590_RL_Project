% June 2010  -  Andrea 
% use the data contained in file .mat to plot convergence in different
% parameters
%
% assume R=1
% data depend essentially on 4 parameters: [ref/h,  xi,  omega,  p]


close all;clear all;
tic
MyFileName = 'PWDG_error_data';
load (MyFileName);

% plotting options
pl1='b-.';  pl2='g--';  pl3='r';
%pl1='bo';  pl2='gd';  pl3='rx';   %use this if PWDG is not ordered





%% collect all the parameters in short vectors
numexp=length(PWDG);
Vp      = [];       % number of equispaced PW directions per element
Vomega  = [];       % wavenumber
VR      = [];       % size of the square domain (R=1)
Vref    = [];       % # of mesh refinements
Vxi     = [];       % parameter xi of the exact solution, integer=regular
Vomegah = [];       % omega * h
Vh      = [];       % meshwidth
Vndof   = [];       % # degrees of freedom, = p * # elements
for i =1:numexp
    Vp = [Vp PWDG{i}.p];
    Vomega  = [Vomega PWDG{i}.omega];
    VR      = [VR PWDG{i}.R];
    Vref    = [Vref PWDG{i}.ref];
    Vxi     = [Vxi PWDG{i}.xi];
    Vomegah = [Vomegah PWDG{i}.omegah];
    Vh      = [Vh PWDG{i}.meshwidth];
    Vndof   = [Vndof PWDG{i}.ndof];
end
Vp      = unique(Vp);
Vomega  = unique(Vomega);
VR      = unique(VR);
Vref    = unique(Vref);
Vxi     = unique(Vxi);
Vomegah = unique(Vomegah);
Vh      = unique(Vh);
Vndof   = unique(Vndof);

%experiments done such that omega and h match and we have few omegah




%% plot analogous to the ones in the p-paper
%omega=8 instead of 10
fix_omega = 8;
fix_ref = 1;

figure
disp('Fig1: convergence wrt p');
xpaper1 =zeros(1,0);ypaper1 =zeros(9,0);
xpaper23=zeros(1,0);ypaper23=zeros(9,0);
xpaper32=zeros(1,0);ypaper32=zeros(9,0);

for i=1:numexp
    if (PWDG{i}.ref ~= fix_ref); continue; end
    if (PWDG{i}.omega ~= fix_omega); continue; end
    if (PWDG{i}.R ~= 1); continue; end
    
    if (PWDG{i}.xi == 1);     
    xpaper1=[xpaper1, PWDG{i}.p];
    ypaper1=[ypaper1, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg 
                        PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg 
                        PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
    %[PWDG{i}.cond_pwdg PWDG{i}.cond_uwvf PWDG{i}.cond_proj]
                    
    elseif (PWDG{i}.xi == 2/3);
    xpaper23=[xpaper23, PWDG{i}.p];
    ypaper23=[ypaper23, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg 
                        PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg 
                        PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];

    elseif (PWDG{i}.xi == 1.5);
    xpaper32=[xpaper32, PWDG{i}.p];
    ypaper32=[ypaper32, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg 
                        PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg 
                        PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
    end
end

subplot(3,3,1)
semilogy(xpaper1,ypaper1(1,:), pl1,'linewidth',2);hold on
semilogy(xpaper1,ypaper1(4,:), pl2,'linewidth',2);
semilogy(xpaper1,ypaper1(7,:), pl3,'linewidth',2);title('xi=1, L2 error');
xlabel ('p');
subplot(3,3,2)
semilogy(xpaper1,ypaper1(2,:), pl1,'linewidth',2);hold on
semilogy(xpaper1,ypaper1(5,:), pl2,'linewidth',2);
semilogy(xpaper1,ypaper1(8,:), pl3,'linewidth',2);title('xi=1, br. H1 error');
subplot(3,3,3)
semilogy(xpaper1,ypaper1(3,:), pl1,'linewidth',2);hold on
semilogy(xpaper1,ypaper1(6,:), pl2,'linewidth',2);
semilogy(xpaper1,ypaper1(9,:), pl3,'linewidth',2);title('xi=1, jumps error');

subplot(3,3,4)
loglog(xpaper23,ypaper23(1,:), pl1,'linewidth',2);hold on
loglog(xpaper23,ypaper23(4,:), pl2,'linewidth',2);
loglog(xpaper23,ypaper23(7,:), pl3,'linewidth',2);title('xi=2/3, L2 error');
subplot(3,3,5)
loglog(xpaper23,ypaper23(2,:), pl1,'linewidth',2);hold on
loglog(xpaper23,ypaper23(5,:), pl2,'linewidth',2);
loglog(xpaper23,ypaper23(8,:), pl3,'linewidth',2);title('xi=2/3, br. H1 error');
subplot(3,3,6)
loglog(xpaper23,ypaper23(3,:), pl1,'linewidth',2);hold on
loglog(xpaper23,ypaper23(6,:), pl2,'linewidth',2);
loglog(xpaper23,ypaper23(9,:), pl3,'linewidth',2);title('xi=2/3, jumps error');

subplot(3,3,7)
loglog(xpaper32,ypaper32(1,:), pl1,'linewidth',2);hold on
loglog(xpaper32,ypaper32(4,:), pl2,'linewidth',2);
loglog(xpaper32,ypaper32(7,:), pl3,'linewidth',2);title('xi=1.5, L2 error');
subplot(3,3,8)
loglog(xpaper32,ypaper32(2,:), pl1,'linewidth',2);hold on
loglog(xpaper32,ypaper32(5,:), pl2,'linewidth',2);
loglog(xpaper32,ypaper32(8,:), pl3,'linewidth',2);title('xi=1.5, br. H1 error');
subplot(3,3,9)
loglog(xpaper32,ypaper32(3,:), pl1,'linewidth',2);hold on
loglog(xpaper32,ypaper32(6,:), pl2,'linewidth',2);
loglog(xpaper32,ypaper32(9,:), pl3,'linewidth',2);title('xi=1.5, jumps error');







%% plot in omegah, to verify pollution, quite useless

select_p = [9, 15];
select_omega = [4, 16, 64];
%select_omega = [2, 8, 32];

figure;
disp('Fig2: convergence wrt omega.h, one line for different p, omega');
for p=select_p
    for omega=select_omega
        xpolh1 =zeros(1,0);ypolh1 =zeros(9,0);
        xpolh23=zeros(1,0);ypolh23=zeros(9,0);
        xpolh32=zeros(1,0);ypolh32=zeros(9,0);
        for i = 1:length(PWDG)
            if (PWDG{i}.p ~= p); continue; end
            if (PWDG{i}.omega ~= omega); continue; end
            
            if (PWDG{i}.xi == 1);
                xpolh1=[xpolh1, PWDG{i}.omegah];
                ypolh1=[ypolh1, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                
            elseif (PWDG{i}.xi == 2/3);
                xpolh23=[xpolh23, PWDG{i}.omegah];
                ypolh23=[ypolh23, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                
            elseif (PWDG{i}.xi == 1.5);
                xpolh32=[xpolh32, PWDG{i}.omegah];
                ypolh32=[ypolh32, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
            end
        end
        
        subplot(3,3,1)
        loglog(xpolh1,ypolh1(1,:), pl1,'linewidth',2);hold on
        loglog(xpolh1,ypolh1(4,:), pl2,'linewidth',2);
        loglog(xpolh1,ypolh1(7,:), pl3,'linewidth',2);title('xi=1, L2 error');
        xlabel ('omegah');
        subplot(3,3,2)
        loglog(xpolh1,ypolh1(2,:), pl1,'linewidth',2);hold on
        loglog(xpolh1,ypolh1(5,:), pl2,'linewidth',2);
        loglog(xpolh1,ypolh1(8,:), pl3,'linewidth',2);title('xi=1, br. H1 error');
        subplot(3,3,3)
        loglog(xpolh1,ypolh1(3,:), pl1,'linewidth',2);hold on
        loglog(xpolh1,ypolh1(6,:), pl2,'linewidth',2);
        loglog(xpolh1,ypolh1(9,:), pl3,'linewidth',2);title('xi=1, jumps error');
                
        subplot(3,3,4)
        loglog(xpolh23,ypolh23(1,:), pl1,'linewidth',2);hold on
        loglog(xpolh23,ypolh23(4,:), pl2,'linewidth',2);
        loglog(xpolh23,ypolh23(7,:), pl3,'linewidth',2);title('xi=2/3, L2 error');
        subplot(3,3,5)
        loglog(xpolh23,ypolh23(2,:), pl1,'linewidth',2);hold on
        loglog(xpolh23,ypolh23(5,:), pl2,'linewidth',2);
        loglog(xpolh23,ypolh23(8,:), pl3,'linewidth',2);title('xi=2/3, br. H1 error');
        subplot(3,3,6)
        loglog(xpolh23,ypolh23(3,:), pl1,'linewidth',2);hold on
        loglog(xpolh23,ypolh23(6,:), pl2,'linewidth',2);
        loglog(xpolh23,ypolh23(9,:), pl3,'linewidth',2);title('xi=2/3, jumps error');
        
        subplot(3,3,7)
        loglog(xpolh32,ypolh32(1,:), pl1,'linewidth',2);hold on
        loglog(xpolh32,ypolh32(4,:), pl2,'linewidth',2);
        loglog(xpolh32,ypolh32(7,:), pl3,'linewidth',2);title('xi=1.5, L2 error');
        subplot(3,3,8)
        loglog(xpolh32,ypolh32(2,:), pl1,'linewidth',2);hold on
        loglog(xpolh32,ypolh32(5,:), pl2,'linewidth',2);
        loglog(xpolh32,ypolh32(8,:), pl3,'linewidth',2);title('xi=1.5, br. H1 error');
        subplot(3,3,9)
        loglog(xpolh32,ypolh32(3,:), pl1,'linewidth',2);hold on
        loglog(xpolh32,ypolh32(6,:), pl2,'linewidth',2);
        loglog(xpolh32,ypolh32(9,:), pl3,'linewidth',2);title('xi=1.5, jumps error');
    end
end



%% plot in omega, again to verify pollution, better

select_p = [15, 25];
select_omegah = Vomegah([3, 6, 9]);

figure;
disp('Fig3: convergence wrt omega, one line for different p, omegah');
for p=select_p
    for omegah=select_omegah
        
        xpol1 =zeros(1,0);ypol1 =zeros(9,0);
        xpol23=zeros(1,0);ypol23=zeros(9,0);
        xpol32=zeros(1,0);ypol32=zeros(9,0);
        for i = 1:length(PWDG)
            if (PWDG{i}.p ~= p); continue; end
            if (PWDG{i}.omegah ~= omegah); continue; end
            
            if (PWDG{i}.xi == 1);
                xpol1=[xpol1, PWDG{i}.omega];
                ypol1=[ypol1, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                
            elseif (PWDG{i}.xi == 2/3);
                xpol23=[xpol23, PWDG{i}.omega];
                ypol23=[ypol23, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                
            elseif (PWDG{i}.xi == 1.5);
                xpol32=[xpol32, PWDG{i}.omega];
                ypol32=[ypol32, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
            end
        end
        
        
        subplot(3,3,1)
        loglog(xpol1,ypol1(1,:), pl1,'linewidth',2);hold on
        loglog(xpol1,ypol1(4,:), pl2,'linewidth',2);
        loglog(xpol1,ypol1(7,:), pl3,'linewidth',2);title('xi=1, L2 error');
        xlabel ('omega');
        subplot(3,3,2)
        loglog(xpol1,ypol1(2,:), pl1,'linewidth',2);hold on
        loglog(xpol1,ypol1(5,:), pl2,'linewidth',2);
        loglog(xpol1,ypol1(8,:), pl3,'linewidth',2);title('xi=1, br. H1 error');
        subplot(3,3,3)
        loglog(xpol1,ypol1(3,:), pl1,'linewidth',2);hold on
        loglog(xpol1,ypol1(6,:), pl2,'linewidth',2);
        loglog(xpol1,ypol1(9,:), pl3,'linewidth',2);title('xi=1, jumps error');
                
        subplot(3,3,4)
        loglog(xpol23,ypol23(1,:), pl1,'linewidth',2);hold on
        loglog(xpol23,ypol23(4,:), pl2,'linewidth',2);
        loglog(xpol23,ypol23(7,:), pl3,'linewidth',2);title('xi=2/3, L2 error');
        subplot(3,3,5)
        loglog(xpol23,ypol23(2,:), pl1,'linewidth',2);hold on
        loglog(xpol23,ypol23(5,:), pl2,'linewidth',2);
        loglog(xpol23,ypol23(8,:), pl3,'linewidth',2);title('xi=2/3, br. H1 error');
        subplot(3,3,6)
        loglog(xpol23,ypol23(3,:), pl1,'linewidth',2);hold on
        loglog(xpol23,ypol23(6,:), pl2,'linewidth',2);
        loglog(xpol23,ypol23(9,:), pl3,'linewidth',2);title('xi=2/3, jumps error');
        
        subplot(3,3,7)
        loglog(xpol32,ypol32(1,:), pl1,'linewidth',2);hold on
        loglog(xpol32,ypol32(4,:), pl2,'linewidth',2);
        loglog(xpol32,ypol32(7,:), pl3,'linewidth',2);title('xi=1.5, L2 error');
        subplot(3,3,8)
        loglog(xpol32,ypol32(2,:), pl1,'linewidth',2);hold on
        loglog(xpol32,ypol32(5,:), pl2,'linewidth',2);
        loglog(xpol32,ypol32(8,:), pl3,'linewidth',2);title('xi=1.5, br. H1 error');
        subplot(3,3,9)
        loglog(xpol32,ypol32(3,:), pl1,'linewidth',2);hold on
        loglog(xpol32,ypol32(6,:), pl2,'linewidth',2);
        loglog(xpol32,ypol32(9,:), pl3,'linewidth',2);title('xi=1.5, jumps error');
    end
end



%% pollution ratio errPWDG/errProj  (= quasi-opt. costant) wrt omega

select_p = [7, 17];     select_omegah = Vomegah([3 5 7 9]);
%select_p = [5:19];      select_omegah = Vomegah([5]);
%select_p = [7];      select_omegah = Vomegah(:);

UseText = 0;  %switch on/off text on the plot

% we could ignore values with higher conditioning, but need to be sure not
% to exclude everything, ... boh, leave it large
cond_threshold = 1e25; 


figure;
disp('Fig4: quasi-opt ratio, omega on the abscissas, one line for different p, omegah');
for p=select_p
    for omegah=select_omegah
        
        xrat1 =zeros(1,0);yrat1 =zeros(3,0);
        xrat23=zeros(1,0);yrat23=zeros(3,0);
        xrat32=zeros(1,0);yrat32=zeros(3,0);
        
        for i = 1:length(PWDG)
            if (PWDG{i}.p ~= p); continue; end
            if (PWDG{i}.omegah ~= omegah); continue; end
            
            if (PWDG{i}.xi == 1);
                if PWDG{i}.cond_pwdg < cond_threshold
                xrat1=[xrat1, PWDG{i}.omega];
                yrat1=[yrat1, [ PWDG{i}.err_pwdg_L2/PWDG{i}.err_proj_L2;  PWDG{i}.err_pwdg_h1/PWDG{i}.err_proj_h1;   PWDG{i}.err_pwdg_dg/PWDG{i}.err_proj_dg ]];
                end
                
            elseif (PWDG{i}.xi == 2/3);
                if PWDG{i}.cond_pwdg < cond_threshold
                xrat23=[xrat23, PWDG{i}.omega];
                yrat23=[yrat23, [ PWDG{i}.err_pwdg_L2/PWDG{i}.err_proj_L2;  PWDG{i}.err_pwdg_h1/PWDG{i}.err_proj_h1;   PWDG{i}.err_pwdg_dg/PWDG{i}.err_proj_dg ]];
                end
                
            elseif (PWDG{i}.xi == 1.5);
                if PWDG{i}.cond_pwdg < cond_threshold                
                xrat32=[xrat32, PWDG{i}.omega];
                yrat32=[yrat32,[ PWDG{i}.err_pwdg_L2/PWDG{i}.err_proj_L2;  PWDG{i}.err_pwdg_h1/PWDG{i}.err_proj_h1;   PWDG{i}.err_pwdg_dg/PWDG{i}.err_proj_dg ]];
                end
            end
        end
        
        
        subplot(3,3,1)
        loglog(xrat1,yrat1(1,:), 'b-o','linewidth',2);hold on;title('xi=1, L2 ratio');
        if (UseText); text(xrat1(1),yrat1(1,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        xlabel ('omega'); ylabel('Pollution Ratio')
        subplot(3,3,2)
        loglog(xrat1,yrat1(2,:), 'b-o','linewidth',2);hold on;title('xi=1, br. H1 ratio');
        if (UseText); text(xrat1(1),yrat1(2,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        subplot(3,3,3)
        loglog(xrat1,yrat1(3,:), 'b-o','linewidth',2);hold on;title('xi=1, jumps ratio');
        if (UseText); text(xrat1(1),yrat1(3,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
               
        subplot(3,3,4)
        loglog(xrat23,yrat23(1,:), 'b-o','linewidth',2);hold on;title('xi=2/3, L2 ratio');
        if (UseText); text(xrat23(1),yrat23(1,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        subplot(3,3,5)
        loglog(xrat23,yrat23(2,:), 'b-o','linewidth',2);hold on;title('xi=2/3, br. H1 ratio');
        if (UseText); text(xrat23(1),yrat23(2,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        subplot(3,3,6)
        loglog(xrat23,yrat23(3,:), 'b-o','linewidth',2);hold on;title('xi=2/3, jumps ratio');
        if (UseText); text(xrat23(1),yrat23(3,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
       
        subplot(3,3,7)
        loglog(xrat32,yrat32(1,:), 'b-o','linewidth',2);hold on;title('xi=1.5, L2 ratio');
        if (UseText); text(xrat32(1),yrat32(1,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        subplot(3,3,8)
        loglog(xrat32,yrat32(2,:), 'b-o','linewidth',2);hold on;title('xi=1.5, br. H1 ratio');
        if (UseText); text(xrat32(1),yrat32(2,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        subplot(3,3,9)
        loglog(xrat32,yrat32(3,:), 'b-o','linewidth',2);hold on;title('xi=1.5, jumps ratio');
        if (UseText); text(xrat32(1),yrat32(3,1), strcat('p=',num2str(p),'oh=',num2str(omegah))); end;
        
    end
end


%% as in the last part of the paper

figure;
disp('Fig5: p-conv for different omega');

select_omega = [0.25 1 4 16 64];
select_ref = 1;

for omega=select_omega
    for ref = select_ref
        
        xpap1 =zeros(1,0);ypap1 =zeros(9,0);
        xpap23=zeros(1,0);ypap23=zeros(9,0);
        xpap32=zeros(1,0);ypap32=zeros(9,0);
        
        for i = 1:length(PWDG)
            if (PWDG{i}.ref ~= ref); continue; end
            if (PWDG{i}.omega ~= omega); continue; end
            
            if (PWDG{i}.xi == 1);
                xpap1=[xpap1, PWDG{i}.p];
                ypap1=[ypap1, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                %[PWDG{i}.cond_pwdg PWDG{i}.cond_uwvf PWDG{i}.cond_proj]
                
            elseif (PWDG{i}.xi == 2/3);
                xpap23=[xpap23, PWDG{i}.p];
                ypap23=[ypap23, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
                
            elseif (PWDG{i}.xi == 1.5);
                xpap32=[xpap32, PWDG{i}.p];
                ypap32=[ypap32, [ PWDG{i}.err_pwdg_L2;  PWDG{i}.err_pwdg_h1;   PWDG{i}.err_pwdg_dg
                    PWDG{i}.err_uwvf_L2;  PWDG{i}.err_uwvf_h1;   PWDG{i}.err_uwvf_dg
                    PWDG{i}.err_proj_L2;  PWDG{i}.err_proj_h1;   PWDG{i}.err_proj_dg  ]];
            end
        end
        
        
        subplot(3,3,1)
        semilogy(xpap1,ypap1(1,:), pl1,'linewidth',2);hold on
        semilogy(xpap1,ypap1(4,:), pl2,'linewidth',2);
        semilogy(xpap1,ypap1(7,:), pl3,'linewidth',2);title('xi=1, L2 error');
        xlabel ('p');
        subplot(3,3,2)
        semilogy(xpap1,ypap1(2,:), pl1,'linewidth',2);hold on
        semilogy(xpap1,ypap1(5,:), pl2,'linewidth',2);
        semilogy(xpap1,ypap1(8,:), pl3,'linewidth',2);title('xi=1, br. H1 error');
        subplot(3,3,3)
        semilogy(xpap1,ypap1(3,:), pl1,'linewidth',2);hold on
        semilogy(xpap1,ypap1(6,:), pl2,'linewidth',2);
        semilogy(xpap1,ypap1(9,:), pl3,'linewidth',2);title('xi=1, jumps error');
        
        subplot(3,3,4)
        loglog(xpap23,ypap23(1,:), pl1,'linewidth',2);hold on
        loglog(xpap23,ypap23(4,:), pl2,'linewidth',2);
        loglog(xpap23,ypap23(7,:), pl3,'linewidth',2);title('xi=2/3, L2 error');
        subplot(3,3,5)
        loglog(xpap23,ypap23(2,:), pl1,'linewidth',2);hold on
        loglog(xpap23,ypap23(5,:), pl2,'linewidth',2);
        loglog(xpap23,ypap23(8,:), pl3,'linewidth',2);title('xi=2/3, br. H1 error');
        subplot(3,3,6)
        loglog(xpap23,ypap23(3,:), pl1,'linewidth',2);hold on
        loglog(xpap23,ypap23(6,:), pl2,'linewidth',2);
        loglog(xpap23,ypap23(9,:), pl3,'linewidth',2);title('xi=2/3, jumps error');
        
        subplot(3,3,7)
        loglog(xpap32,ypap32(1,:), pl1,'linewidth',2);hold on
        loglog(xpap32,ypap32(4,:), pl2,'linewidth',2);
        loglog(xpap32,ypap32(7,:), pl3,'linewidth',2);title('xi=1.5, L2 error');
        subplot(3,3,8)
        loglog(xpap32,ypap32(2,:), pl1,'linewidth',2);hold on
        loglog(xpap32,ypap32(5,:), pl2,'linewidth',2);
        loglog(xpap32,ypap32(8,:), pl3,'linewidth',2);title('xi=1.5, br. H1 error');
        subplot(3,3,9)
        loglog(xpap32,ypap32(3,:), pl1,'linewidth',2);hold on
        loglog(xpap32,ypap32(6,:), pl2,'linewidth',2);
        loglog(xpap32,ypap32(9,:), pl3,'linewidth',2);title('xi=1.5, jumps error');
        
     end
end



%% Ilaria:
select_p = Vp(3);       select_h=  Vh(1);       select_omega=Vomega(1:10);
%select_p = Vp(1);       select_h=  Vh(1);       select_omega=Vomega(1:7);

% we could ignore values with higher conditioning, but need to be sure not
% to exclude everything, ... boh, leave it large
cond_threshold = 1e25; 
p=select_p;
h=select_h;
xrat1=[];yrat1=[];yerr1=[];ybestappr1=[];
for omega=select_omega
 for i = 1:length(PWDG)
   if (PWDG{i}.p ~= p); continue; end
   if (PWDG{i}.omega ~= omega); continue; end
   if (PWDG{i}.xi ~= 1); continue; end
   if (PWDG{i}.meshwidth ~= h); continue; end
   
   if PWDG{i}.cond_pwdg < cond_threshold 
         xrat1=[xrat1, PWDG{i}.omega];
         yrat1=[yrat1, PWDG{i}.err_pwdg_L2/PWDG{i}.err_proj_L2];
         yerr1=[yerr1, PWDG{i}.err_pwdg_L2];
         ybestappr1=[ybestappr1, PWDG{i}.err_proj_L2];         
   end
 end

end
format long e
select_omega
xrat1
yrat1
yerr1
ybestappr1

figure;
loglog(xrat1,yrat1(1,:), 'b-o','linewidth',2);hold on;title('xi=1, L2 ratio');
xlabel ('omega'); ylabel('Pollution Ratio')



toc
