


clc;clear all;close all;

scriptLocation = fileparts(fileparts(mfilename('fullpath') ));
addpath([scriptLocation filesep 'scripts']);
addpathFolderStructure()

load(['results' filesep 'DataMatTot_MacPcCombined'])
load(['data' filesep 'elXY'])

par.rmodes = 25;
%%
sensor_vec = [];
sensorMatTotSq = squeeze(sensorMatTot(2,:,:,:));
index_mat = [];
for j = 1:length(sensorMatTotSq)
    n_iter = length(nonzeros(sensorMatTotSq(j,1,:)));
    for k = 1:n_iter 
        pre_length = length(sensor_vec);
        add_length = length( sensorMatTotSq(j,1:j,k));
        sensor_vec = [sensor_vec , sensorMatTotSq(j,1:j,k) ];
        index_mat(j,k,1:add_length)= pre_length + (1:add_length) ;
    end
end

sensor_x = x(sensor_vec);
sensor_y = y(sensor_vec);

%%
% rad_list = [0.001,0.01,0.1];
rad_list = [0.0005,0.005,0.05,0.5];
% rad_list = [1,10];
col = linspecer(length(rad_list));
tic 
%     figure();
% [xG,yG] = meshgrid( -1.25:0.0005:1.25,0:0.0005:5);
[xG,yG] = meshgrid( -1.25:0.0005:1.25,0:0.0005:5);
rad_mat = [];
scatter_vec = [];



% C = cumsum(gauss(:))/sum(gauss(:));
% placeInd = 1+sum(C(end)*rand>C);


    figure()
for jj= 1:length(rad_list) 
    rad = rad_list(jj)
%     figure()
    count = 1;
    for j = 1:length(sensor_vec)
        
        gauss = exp( - ((xG-sensor_x(j)).^2 + (yG-sensor_y(j)).^2)/rad_list(jj));
        
        X = xG(:);
        Y = yG(:);
        C = cumsum(gauss(:))/sum(gauss(:));
        placeInd = 1+sum(C(end)*rand>C);
        x_place(count) = X(placeInd );
        y_place(count) = Y(placeInd );
        count = count + 1;
% 
        scatter(sensor_x(j),sensor_y(j),'x')
        hold on
        b = scatter(x_place(j),y_place(j),10,col(jj,:),'fill');
    end
    scatter_vec = [scatter_vec ; b];
    rad_mat = [ rad_mat , jj*ones(1,length(sensor_vec)) ];
end
loc_error = [x_place;y_place];
toc 
axis equal
legend(scatter_vec,num2str(rad_list'))

% figure();
% scatter(sensor_x(1:j),sensor_y(1:j))
% hold on
% scatter(x_place(1:j),y_place(1:j),'.')

%%
strainSet = eulerLagrangeConcatenateXYloc( 0,0, sensor_vec,loc_error, par); 

% apply neural filter
[X,G] = neuralEncoding(strainSet, par );

%% 
for jj= 1:length(rad_list) 
%     X(rad_mat == jj,:)
    [Xtrain, Xtest, Gtrain, Gtest] = rand_cross_val(X(rad_mat == jj,:), G, par.trainFraction);

        for j = 1:length(sensorMatTotSq)
            n_iter = length(nonzeros(sensorMatTotSq(j,1,:)));
            for k = 1:n_iter 
                accuracy = sensorLocClassify(  nonzeros(index_mat(j,k,:)) ,Xtrain,Gtrain,Xtest,Gtest );

                DataMat(jj,j,k) =accuracy; 
            end
        end
end



save(['results' filesep 'accuracy_after_error.mat'],'DataMat')


