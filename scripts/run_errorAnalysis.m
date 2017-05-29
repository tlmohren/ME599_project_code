


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
rad_list = [0.005,0.05,0.5];
tic 
%     figure();
[xG,yG] = meshgrid( -1.25:0.005:1.25,0:0.005:5);
count = 1;
rad_mat = [];

for jj= 1:length(rad_list) 
    rad = rad_list(jj);
    for j = 1:length(sensor_vec)

        gauss = exp( - ((xG-sensor_x(j)).^2 + (yG-sensor_y(j)).^2)/rad);
        C = cumsum(gauss(:))/sum(gauss(:));
        placeInd = 1+sum(C(end)*rand>C);
        x_place(count) = xG(placeInd );
        y_place(count) = yG(placeInd );
        count = count + 1;
    %     surf(xG,yG,gauss)
    %     hold on
    %     scatter3(x_place(j),y_place(j),1,'fill','r')
    end
    rad_mat = [ rad_mat , jj*ones(1,length(sensor_vec)) ];
end
loc_error = [x_place;y_place];
toc 

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


