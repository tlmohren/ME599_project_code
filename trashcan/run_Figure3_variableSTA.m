
%% Find nonzero elements in matrix 
clc;clear all; close all 
scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
addpath([scriptLocation filesep 'scripts']);
addpathFolderStructure()

% load(['results' filesep 'analysis_FigR1toR4_XXYY_270Par'])

load(['results' filesep 'tempDataMatTot'])
Datamat = dataMatTot;

%% 

% Datamat = rand(size(Datamat))*0.2 + 0.6;

%% see which simulations belong to this parameter set 
test = ( [varParList.phi_dist] == 0)& ...
            ( [varParList.theta_dist] == 0);
figure();plot(test)


%% 
% bin_SSPOCon = ( [varParList.phi_dist] == 0) & ...
%             ( [varParList.theta_dist] == 0) & ...
%             ( [varParList.SSPOCon] == 1 ) & ...
%             ( [varParList.xInclude] == 1) & ...%%%%% note error in assignment of yinlude~!!
%             ( [varParList.yInclude] == 1) & ...%%%%% note error in assignment ofyinlude~!!
%             ( [varParList.NLDshift] == 0.5) & ...
%             ( [varParList.NLDsharpness] == 10);
% bin_SSPOCoff = ( [varParList.phi_dist] == 0)  & ...
%             ( [varParList.theta_dist] == 0) & ...
%             ( [varParList.SSPOCon] == 0 ) & ...
%             ( [varParList.xInclude] == 1) & ...    %%%%% note error in assignment of yinlude~!!
%             ( [varParList.yInclude] == 1) & ...%%%%% note error in assignment ofyxinlude~!!
%             ( [varParList.NLDshift] == 0.5) & ...
%             ( [varParList.NLDsharpness] == 10);
%         
        
        
bin_SSPOCon = ( [varParList_short.phi_dist] == 0) & ...
            ( [varParList_short.theta_dist] == 0) & ...
            ( [varParList_short.SSPOCon] == 1 ) & ...
            ( [varParList_short.xInclude] == 0) & ...%%%%% note error in assignment of yinlude~!!
            ( [varParList_short.yInclude] == 1) & ...%%%%% note error in assignment ofyinlude~!!
            ( [varParList_short.NLDshift] == 0.5) & ...
            ( [varParList_short.NLDsharpness] == 10);
bin_SSPOCoff = ( [varParList_short.phi_dist] == 0)  & ...
            ( [varParList_short.theta_dist] == 0) & ...
            ( [varParList_short.SSPOCon] == 0 ) & ...
            ( [varParList_short.xInclude] == 0) & ...    %%%%% note error in assignment of yinlude~!!
            ( [varParList_short.yInclude] == 1) & ...%%%%% note error in assignment ofyxinlude~!!
            ( [varParList_short.NLDshift] == 0.5) & ...
            ( [varParList_short.NLDsharpness] == 10);
        
        
        
        
        
        
        
        
        
        
        
figure();plot(bin_SSPOCon,'-o')
hold on;plot(bin_SSPOCoff,'-o')

ind_SSPOCon = find(bin_SSPOCon);
ind_SSPOCoff = find(bin_SSPOCoff);

par.STAwidthList = [1:2:10];
par.STAshiftList = [-11:-2:-20];% 


        par.STAwidthList = [1:2:10];
        par.STAshiftList = [-4:-3:-16];% 
        par.NLDshiftList = [-0.1:0.2:0.9];
        par.NLDsharpnessList = [4:2:12];
        
        
        
% n_plots = 16; 
n_x = (length(par.STAwidthList) + 1);
n_y = (length(par.STAshiftList) +1); 
n_plots = n_x*n_y; 
% col = fliplr({'-k','-r'});
% dotcol = fliplr({'.k','.r'}); 
% col = ({'-k','-r'});
col = {ones(3,1)*0.5,'-r'};
dotcol = ({'.k','.r'}); 

%% create subplot routine


mat_1 = reshape(1:n_x*n_y,[n_x,n_y])';
mat_2 = mat_1(2:end,2:end);
vec_3 = sort(mat_2(:));
vec_v = mat_1(2:end,1);

% color_vec = {'-r','-k'};
fig2=figure('Position', [100, 100, 1000, 800]);
for j = 1:length(vec_3)
    subplot(n_y,n_x, vec_3(j))
%     j 
%     vec_3(j)
    % plot sspoc on 
        Dat_I = ind_SSPOCoff(j);
        for k2 = 1:size(Datamat,2)
            meanVec(k2) = mean(  nonzeros(Datamat(Dat_I,k2,:))   );
            stdVec(k2) = std(  nonzeros(Datamat(Dat_I,k2,:))   );
% 
            iters = length(nonzeros(Datamat(Dat_I,k2,:)) );
%             scatter( ones(iters,1)*k2,nonzeros(Datamat(Dat_I,k2,:)) , dotcol{1})
            hold on
   
        end
        realNumbers = find(~isnan(meanVec));
        a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{1},0.8);
    
   % plot sspoc om
        Dat_I = ind_SSPOCon( j);
        for k2 = 1:size(Datamat,2)
            iters = length(nonzeros(Datamat(Dat_I,k2,:)) );
            if iters > 0
                meanVec(k2) = mean(  nonzeros(Datamat(Dat_I,k2,:))   );
                stdVec(k2) = std(  nonzeros(Datamat(Dat_I,k2,:))   );
            else
                meanVec(k2) = nan;
                stdVec(k2) = nan;
            end
            scatter( ones(iters,1)*k2,nonzeros(Datamat(Dat_I,k2,:)) , dotcol{2})
            hold on
        end
        realNumbers = find(~isnan(meanVec));
%         a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{2},0.8);
        
        plot(realNumbers,meanVec(realNumbers),col{2})
        
%         yticks([0.4:0.2:1])
        ax = gca;
        ax.YTick = 0.4:0.2:1;
        axis([0,30,0.4,1])
%         axis off
end

for j = 2:n_x
   subplot(n_y,n_x,j )
   par.STAwidth = par.STAwidthList(j-1);
   par.STAshift = -10;
   t_sta = -39:0;
      par.STAFunc = @(t)  2 * exp( -(t-par.STAshift) .^2 ...
            ./ (2*par.STAwidth ^2) ) ...
            ./ (sqrt(3*par.STAwidth) *pi^1/4)...
            .* ( 1-(t-par.STAshift).^2/par.STAwidth^2);
        par.STAfilt = par.STAFunc(t_sta);   
        

            if par.STAwidth == 3
                plot(par.STAfilt,'k','LineWidth',4); hold on;
            else
                plot(par.STAfilt); hold on;
            end
   axis off
end

%% 
for j = 1:length(vec_v)
   subplot(n_y,n_x,vec_v(j) )
%    par.STAwidth = par.STAwidthList(1);
   par.STAshift = par.STAshiftList(j);

    for k = 1:length(par.STAwidthList)
        par.STAwidth = par.STAwidthList(k);
        
   t_sta = -39:0;
      par.STAFunc = @(t)  2 * exp( -(t-par.STAshift) .^2 ...
            ./ (2*par.STAwidth ^2) ) ...
            ./ (sqrt(3*par.STAwidth) *pi^1/4)...
            .* ( 1-(t-par.STAshift).^2/par.STAwidth^2);
        par.STAfilt = par.STAFunc(t_sta);   
%         plot(t_sta,par.STAfilt/max(par.STAfilt),'Color',ones(1,3)*0.7,'LineWidth',0.5);hold on
        
            if par.STAshift == -10 && par.STAwidth == 3
                plot(par.STAfilt,'k','LineWidth',4); hold on;
             plot(t_sta,par.STAfilt/max(par.STAfilt),'k','LineWidth',2);hold on
            else
             plot(t_sta,par.STAfilt/max(par.STAfilt),'Color',ones(1,3)*0.7,'LineWidth',0.5);hold on
            end
            
            
            
    end
   x_mid = mean(par.STAshiftList);
   plot([x_mid,x_mid],[-0.5,1.2],'Color',ones(1,3)*0.5,'LineWidth',[2])
   hold on 

   quiver( x_mid,0.5,  - x_mid + par.STAshiftList(j),0,'m','LineWidth',3,'MaxHeadSize',1)


   plot([0,0]+ par.STAshiftList(j),[-0.5,1.2] ,'k','LineWidth',[2])

   axis off
   axis([-39,0,-0.5,1.2])
end
%% 
if 0 
    saveas(fig2,['figs' filesep 'Figure3_variableSTA'], 'png')
end