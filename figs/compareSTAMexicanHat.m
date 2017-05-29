%% 

clc;clear all;close all

%% 

        par.STAwidthList = [1:2:10];
        par.STAshiftList = [-4:-3:-16];% 
        par.NLDshiftList = [-0.1:0.2:0.9];
        par.NLDsharpnessList = [4:2:12];

fig1 = figure('Position', [100, 100, 1200, 600]);


subplot(231)
load('exp_STA');
plot(sta_1kHz/max(sta_1kHz));


    par.STAwidth = 3;% par.STAwidthList(j-1);
   par.STAshift = -10;%par.STAshiftList(1);
   t_sta = -39:0;
      par.STAFunc = @(t)  2 * exp( -(t-par.STAshift) .^2 ...
            ./ (2*par.STAwidth ^2) ) ...
            ./ (sqrt(3*par.STAwidth) *pi^1/4)...
            .* ( 1-(t-par.STAshift).^2/par.STAwidth^2);
        par.STAfilt = par.STAFunc(t_sta);   
 hold on;
plot(par.STAfilt)
legend('Experimental STA','Mexican hat','Location','Best')
% legend('Experimental STA','Mexican hat','Location','NorthEastOutside')

subplot(232)
%         par.STAwidthList = [1:2:10];
        par.STAwidth = 3;
        for j= 1:length(par.STAwidthList)
            par.STAwidth = par.STAwidthList(j);
            par.STAFunc = @(t)  2 * exp( -(t-par.STAshift) .^2 ...
                ./ (2*par.STAwidth ^2) ) ...
                ./ (sqrt(3*par.STAwidth) *pi^1/4)...
                .* ( 1-(t-par.STAshift).^2/par.STAwidth^2);
            par.STAfilt = par.STAFunc(t_sta);   
            if par.STAwidth == 3
                plot(par.STAfilt,'k','LineWidth',4); hold on;
            else
                plot(par.STAfilt,'Color',ones(1,3)*0.7); hold on;
            end
        end
            
%         par.STAshiftList = [-1:-2:-10];% +
subplot(233)
%         par.STAshiftList = [-4:-3:-16];% 
        par.STAwidth = 3;
        for j= 1:length(par.STAshiftList)
            par.STAshift = par.STAshiftList(j);
            par.STAFunc = @(t)  2 * exp( -(t-par.STAshift) .^2 ...
                ./ (2*par.STAwidth ^2) ) ...
                ./ (sqrt(3*par.STAwidth) *pi^1/4)...
                .* ( 1-(t-par.STAshift).^2/par.STAwidth^2);
            par.STAfilt = par.STAFunc(t_sta);   
            if par.STAshift == -10
                display('happened')
                plot(par.STAfilt,'k','LineWidth',4); hold on;
            else
                plot(par.STAfilt,'Color',ones(1,3)*0.7); hold on;
            end
        end
%         
%         
%%  
% 
load('exp_NLD')
subplot(234) 
    par.NLDshift = 0.4813;
    par.NLDsharpness = 8.42;
    % par.NLDsharpness = par.NLDsharpnessList(j-1);
        par.NLD = @(s) 1./(  1 +...
            exp( -(s-par.NLDshift) * par.NLDsharpness)  );
        x = -1:0.02:1;
    plot(x,par.NLD(x))
    legend('Experimental NLD','Location','Best')
    % legend('Experimental NLD','Location','NorthEastOutside')
        
        
     
%         par.NLDshiftList = [-0.2:0.2:0.8];
%         par.NLDsharpnessList = [5:2:15];
    
subplot(235) 
    par.NLDshift = 0.4813;
    
    for j = 1:length(par.NLDsharpnessList)
        par.NLDsharpness = par.NLDsharpnessList(j);
        % par.NLDsharpness = par.NLDsharpnessList(j-1);
            par.NLD = @(s) 1./(  1 +...
                exp( -(s-par.NLDshift) * par.NLDsharpness)  );
            x = -1:0.02:1;
            if  par.NLDsharpness  == 8
%                 display('happened')
                plot(x,par.NLD(x),'k','LineWidth',4); hold on;
            else
                plot(x,par.NLD(x),'Color',ones(1,3)*0.7); hold on;
            end
    end
%         legend('Experimental NLD','Location','Best')
    % legend('Experimental NLD','Location','NorthEastOutside')
    
subplot(236) 
    par.NLDsharpness = 8.42;
    
    for j = 1:length(par.NLDshiftList)
        par.NLDshift = par.NLDshiftList(j);
        % par.NLDsharpness = par.NLDsharpnessList(j-1);
            par.NLD = @(s) 1./(  1 +...
                exp( -(s-par.NLDshift) * par.NLDsharpness)  );
            x = -1:0.02:1;
            if  par.NLDshift  == 0.5
%                 display('happened')
                plot(x,par.NLD(x),'k','LineWidth',4); hold on;
            else
                plot(x,par.NLD(x),'Color',ones(1,3)*0.7); hold on;
            end
    end
    
    
    %%
    
    
saveas(fig1,['Figure03_filterVariation'], 'png')
    
    