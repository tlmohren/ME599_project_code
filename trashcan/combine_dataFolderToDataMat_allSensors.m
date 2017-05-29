% analyze varParList 
clc;clear all;close all
%% 
load(  ['data',filesep, 'ParameterList_allSensors']  )
par.varParNames = fieldnames(varParList);
par.rmodes = 30;

%% 
dataMatTot = zeros( length(varParList) ,par.iter);


saveNameCell = {};
saveNameCount = 0;
duplicates = 0;
%% 
aa = dir(['data' filesep  'Data*']);
% for j = 1
counter= 0
for j = 1:length(varParList)
    for k = 1:length(par.varParNames)
        par.(par.varParNames{k}) = varParList(j).(par.varParNames{k});
    end
    
%         length(nonzeros(dataMatTot))
    for j2 = 1326
%      
        saveName = sprintf('Data_dT%g_dP%g_xIn%g_yIn%g_sOn%g_STAw%g_STAs%g_NLDs%g_NLDg%g_wT%g',...
                            [par.theta_dist , par.phi_dist , par.xInclude , par.yInclude , par.SSPOCon , ...
                            par.STAwidth , par.STAshift , par.NLDshift , par.NLDsharpness , j2 ]);
%          if %%
             duplic_check = 0;
             for j4 = 1:length(saveNameCell)
                 if length(saveNameCell{j4}) == length(saveName) 
                     if   saveNameCell{j4} == saveName; 
                         duplic_check = duplic_check + 1;
                     end
                 end
             end
             saveNameCell{saveNameCount + 1} = saveName;
             saveNameCount = saveNameCount+1;
%          else
%              duplicates = duplicates + 1;
%          end
          bb = dir(['data' filesep saveName '_*']);
          
          if length(bb)>0
              for k2 = 1:length(bb)
                 
                a3 = load( ['data' filesep bb(k2).name] ); 
                
                
                
                [I,J] = ind2sub(size(a3.DataMat),find(a3.DataMat));
                  
                  for j3 = 1:length(I)
                        q = I(j3);
            %             a3.DataMat(I(j3),J(j3))
                        prev = length(find( dataMatTot(j,:) )  );
                        dataMatTot(j, prev+1) = a3.DataMat(I(j3),J(j3)); 
                  end
              end
          else
                counter= counter+1;
          end
      
    end
      
      
end

%% save datamattot

save(['results' filesep  'tempDataMatTot_allSensors.mat'],'dataMatTot','varParList','varParList','par')
 