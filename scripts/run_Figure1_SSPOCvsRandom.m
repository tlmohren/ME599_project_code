% generate R1 figure, just 

clc;clear all;close all;

scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
addpath([scriptLocation filesep 'scripts']);
addpathFolderStructure()
%% 

load(['results' filesep 'DataMatTot_MacPcCombined'])
jumbled = load(['results' filesep 'accuracy_after_error.mat']);
noisy_Datamat = jumbled.DataMat;
rad_list = jumbled.rad_list;
Datamat = dataMatTot;

legend_vec = [];
col = {'-k','-r','-b','-m','-y'};
dotcol = {'.k','.r','.b','.m','.y'}; 
fig1 = figure('Position', [100, 100, 1000, 800]);
%% random 
j = 1;
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
col = linspecer(size(noisy_Datamat,1)+1);
j = 2;
clear meanVec stdVec
for k = 1:size(Datamat,2)
    meanVec(k) = mean(  nonzeros(Datamat(j,k,:))   );
    stdVec(k) = std(  nonzeros(Datamat(j,k,:))   );

    iters = length(nonzeros(Datamat(j,k,:)) );
%         scatter( ones(iters,1)*k,nonzeros(Datamat(j,k,:)) , dotcol{j})
    hold on

end
realNumbers = find(~isnan(meanVec)); 
b = plot(realNumbers,meanVec(realNumbers),'Color',col(1,:));
%     a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{j},0.8);

legend_vec = [legend_vec,b];
%% Error 

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
    b = plot( realNumbers,meanVec(realNumbers), 'Color',col(j+1,:));
%         a = shadedErrorBar(realNumbers, meanVec(realNumbers),stdVec(realNumbers),col{j},0.8);
legend_vec = [legend_vec,b];
end

%%
names = {'Random sensor placement', 'Optimal placement'};
for j = 1:length(rad_list)
   names{j+2} = ['\sigma = ',num2str(rad_list(j)/2.5*100) ,'% of chord']; 
end
legend(legend_vec,names,'Location','SouthEast')

axis([0,30,0.4,1])
xlabel('\# sensors')
ylabel('Accuracy [-]')
grid on

saveas(fig1,['figs' filesep 'Figure1_SSPOCvsRandom_withlocation_error'], 'png')

