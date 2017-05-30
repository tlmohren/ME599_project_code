% code for error analysis 
clc;clear all;close all;

scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
addpath([scriptLocation filesep 'scripts']);
addpathFolderStructure()
load(['results' filesep 'DataMatTot_MacPcCombined'])
load(['data' filesep 'elXY'])

x = x';
y = y';

%%
sensor_vec = [];
sensorMatTotSq = squeeze(sensorMatTot(2,:,:,:));
index_conversion = [];
for j = 1:size(sensorMatTotSq,1) % what w_trunc 
    n_iter = length(nonzeros(sensorMatTotSq(j,1,:)));
    for k = 1:n_iter 
        pre_length = length(sensor_vec);
        add_length = length( sensorMatTotSq(j,1:j,k));
        
        sensors_added = sensorMatTotSq(j,1:j,k);
        sensors_index = [ones(1,length(sensors_added))*j; ones(1,length(sensors_added))*k];
        
        if nonzeros(sensors_added) >0
            sensor_vec = [sensor_vec , sensors_added];
            index_conversion = [index_conversion , sensors_index ];
        end
    end
end

sensor_x = x(sensor_vec);
sensor_y = y(sensor_vec);

%%
rad_list = [0.1,0.5,1,2.5];
col = linspecer(length(rad_list));
[xG,yG] = meshgrid( -1.25:0.01:1.25,0:0.01:5);
rad_mat = [];
scatter_vec = [];

tic 
count = 1;
        
figure()
for jj= 1:length(rad_list) 
    rad = rad_list(jj);
    for j = 1:300%length(sensor_vec)
        gauss = exp( - 0.5*((xG-sensor_x(j)).^2 + (yG-sensor_y(j)).^2)/ rad_list(jj)^2);
%         
        C = cumsum(gauss(:))/sum(gauss(:));
        placeInd = 1+sum(C(end)*rand>C);
        x_place(count) = xG(placeInd );
        y_place(count) = yG(placeInd );
% 
         scatter(sensor_x(j),sensor_y(j),'xk');
        hold on
        b = scatter(x_place(count),y_place(count),10,col(jj,:),'fill');
        
        count = count + 1;
    end
    scatter_vec = [scatter_vec ; b];
    rad_mat = [ rad_mat , jj*ones(1,length(sensor_vec)) ];
end
toc 

a = scatter(sensor_x(j),sensor_y(j),'xk');
scatter_vec = [scatter_vec ; a];
loc_error = [x_place;y_place];

for j = 1:length(rad_list)
    legend_names{j} = ['\sigma = ',num2str(rad_list(j)/2.5*100) ,'% of chord'];
end
legend_names{length(rad_list) +1} = 'optimal sensor locations';
axis equal
legend(scatter_vec,legend_names,'Location','NorthEastOutside')

%%
% figure();
%         gauss = exp( - 0.5*((xG-sensor_x(j)).^2 + (yG-sensor_y(j)).^2)/ 0.25^2);
% surf(xG,yG,gauss)
% axis equal
% scatter(sensor_x(1:j),sensor_y(1:j))
% hold on
% scatter(x_place(1:j),y_place(1:j),'.')

%% 
strainSet = eulerLagrangeConcatenateXYloc( 0,0, sensor_vec,loc_error, par); 
% apply neural filter
[X,G] = neuralEncoding(strainSet, par );

%% 
for jj= 1:length(rad_list) 
    [Xtrain, Xtest, Gtrain, Gtest] = rand_cross_val(X(rad_mat == jj,:), G, par.trainFraction);
    for j = 1:length(sensorMatTotSq)
        n_iter = length(nonzeros(sensorMatTotSq(j,1,:)));
        for k = 1:n_iter 
            sensors_noise = sensor_vec( (index_conversion(1,:) == j) &...
            (index_conversion(2,:) == k) );
            accuracy = sensorLocClassify(  sensors_noise,Xtrain,Gtrain,Xtest,Gtest );
            DataMat(jj,j,k) =accuracy; 
        end
    end
end

%% 
save(['results' filesep 'accuracy_after_error.mat'],'DataMat','rad_list')


