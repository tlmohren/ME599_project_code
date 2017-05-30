% generate R1 figure, just 

clc;clear all;close all;

scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
addpath([scriptLocation filesep 'scripts']);
addpathFolderStructure()
%% 

% load(['results' filesep 'analysis_FigR1toR4_yOnly_87Par.mat'])


% load(['results' filesep 'analysis_FigR1toR4_XXYY_270Par'])

% load(['results' filesep 'tempDataMatTot'])

load(['results' filesep 'DataMatTot_MacPcCombined'])
jumbled = load(['results' filesep 'accuracy_after_error.mat']);
noisy_Datamat = jumbled.DataMat;

Datamat = dataMatTot;



legend_vec = [];



col = {'-k','-r','-b','-m','-y'};
dotcol = {'.k','.r','.b','.m','.y'}; 
fig1 = figure('Position', [100, 100, 1000, 800]);
% for j = 1:2
%% random 
    j = 1
            clear meanVec stdVec
    for k = 1:size(Datamat,2)
        meanVec(k) = mean(  nonzeros(Datamat(j,k,:))   );
        stdVec(k) = std(  nonzeros(Datamat(j,k,:))   );

        iters = length(nonzeros(Datamat(j,k,:)) );
%         scatter( ones(iters,1)*k,nonzeros(Datamat(j,k,:)) , dotcol{j})
        hold on
    
    end
    realNumbers = find(~isnan(meanVec));
    a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{j},0.8);
    legend_vec = [legend_vec,a.mainLine];
%% SSPOC
    j = 2
            clear meanVec stdVec
    for k = 1:size(Datamat,2)
        meanVec(k) = mean(  nonzeros(Datamat(j,k,:))   );
        stdVec(k) = std(  nonzeros(Datamat(j,k,:))   );

        iters = length(nonzeros(Datamat(j,k,:)) );
%         scatter( ones(iters,1)*k,nonzeros(Datamat(j,k,:)) , dotcol{j})
        hold on
    
    end
    realNumbers = find(~isnan(meanVec)); 
    b = plot(realNumbers,meanVec(realNumbers), col{j});
    a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{j},0.8);

    legend_vec = [legend_vec,b];
%% Error 
dotcol = {'.b','.m','.g','.y','.m','.g'};  
col = {'-b','-m','-g','-y','-m','-g'}; 
%     j = 3
    for j = 1:size(noisy_Datamat,1)
%         figure()
            clear meanVec stdVec
        for k = 1:size(noisy_Datamat ,2)
            meanVec(k) = mean(  nonzeros(noisy_Datamat(j,k,:))   );
            stdVec(k) = std(  nonzeros(noisy_Datamat(j,k,:))   );

            iters = length(nonzeros(noisy_Datamat(j,k,:)) );
%             scatter( ones(iters,1)*k,nonzeros(noisy_Datamat(j,k,:)) , dotcol{j})
            hold on
        end
        realNumbers = find(~isnan(meanVec)); 
%         figure()
        b = plot( realNumbers,meanVec(realNumbers), col{j});
        a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{j},0.8);
    legend_vec = [legend_vec,b];
    end

       
    
    
    legend(legend_vec,{'Random sensor placement', 'Optimal placement', 'placement with 0.001 % error', 'placement with 0.01 % error',...
                'placement with 0.1 % error', 'placement with 1 % error'},'Location','SouthEast')
    
    
    
    
    
% end
axis([0,30,0.4,1])
xlabel('\# sensors')
ylabel('Accuracy [-]')
grid on

saveas(fig1,['figs' filesep 'Figure1_SSPOCvsRandom_withlocation_error'], 'png')

%% see which simulations belong to this parameter set 

% 
% varParCase = 2;
% q_select = 16;
%         n_iters= length(nonzeros(Datamat(varParCase,q_select,:)))
%         
% %         length(nonzeros(sensorMatTot(varParCase,q_select,:,:)))
%         binar = zeros(26*51,1);
%         for j = 1:n_iters
% %             sensorMatTot(varParCase,q_select,:,j)
%             binar(sensorMatTot(varParCase,q_select,1:q_select,j)) = binar(sensorMatTot(varParCase,q_select,1:q_select,j)) +1;
%         end
%         binar = binar/n_iters;
%         figure()
%         plotSensorLocs(binar,par)
% %         meanVec(k) = mean(  nonzeros(Datamat(j,k,:))   );
% %         stdVec(k) = std(  nonzeros(Datamat(j,k,:))   );
% % 
% %         iters = length(nonzeros(Datamat(j,k,:)) );
% % %         scatter( ones(iters,1)*k,nonzeros(Datamat(j,k,:)) , dotcol{j})
% %         hold on
% %     
% %     end









