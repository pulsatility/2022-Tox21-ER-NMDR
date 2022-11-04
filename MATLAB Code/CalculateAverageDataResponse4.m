%% Calculate the Average Response and Viability for Chemicals
function [averaged_response_data_with_SID_exclude_concentration,responseDataForReplicatedExperimentsWithPubIDRankOrder,...
    averaged_viability_data_with_SID_exclude_concentration,viabilityDataForReplicatedExperimentsWithPubIDRankOrder,...
    averagedConcentrationDataForResponseWithSID,averagedConcentrationDataForViabilityWithSID] = CalculateAverageDataResponse4(Data_cells)

% response data starts from column 15 to column 29, a total of 15 columns.
% concentration data starts from column 31 to column 45
% pubchem SID: column 49
pubchemSIDColumn = Data_cells{49};

sampleDataType = string(Data_cells{4});
 
rankColumn = [1:1:length(pubchemSIDColumn)]';
 
responseDataColumns = [Data_cells{15:29}];
 
concentrationDataColumns = [Data_cells{31:45}];

responseDataWithRankAndPubchemSID = [rankColumn,pubchemSIDColumn,responseDataColumns,concentrationDataColumns];
% sort matrix based on the second column:pubchemSID
responseDataWithRankAndPubchemSID_sorted = sortrows(responseDataWithRankAndPubchemSID,2);
 
responseDataForReplicatedExperimentsWithPubIDRankOrder = [];
viabilityDataForReplicatedExperimentsWithPubIDRankOrder = [];

responseDataSetForSingleChemical = [];
concentrationForResponseForSingleChemical = [];

viabilityDataSetForSingleChemical = [];
concentrationForViabilityForSingleChemical = [];

averaged_response_data_with_SID_exclude_concentration = [];
averagedResponseData = [];
averagedConcentrationDataForResponse = [];
averagedConcentrationDataForResponseWithSID = [];

averaged_viability_data_with_SID_exclude_concentration = [];
averagedViabilityData = [];
averagedConcentrationDataForViability = [];
averagedConcentrationDataForViabilityWithSID = [];

j = 1;
 
tic
for i = 1:1:length(responseDataWithRankAndPubchemSID_sorted(:,1))
    if (i > 1)
        if (responseDataWithRankAndPubchemSID_sorted(i,2) == responseDataWithRankAndPubchemSID_sorted(i-1,2))
            
            if (contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "antagonist"))
                responseDataForReplicatedExperimentsWithPubIDRankOrder = [responseDataForReplicatedExperimentsWithPubIDRankOrder;...
                    responseDataWithRankAndPubchemSID_sorted(i,:)];
                responseDataSetForSingleChemical = [responseDataSetForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,3:17)];
                concentrationForResponseForSingleChemical = [concentrationForResponseForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,18:32)];
            end
            
            if(contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "viability"))
                viabilityDataForReplicatedExperimentsWithPubIDRankOrder = [viabilityDataForReplicatedExperimentsWithPubIDRankOrder;...
                    responseDataWithRankAndPubchemSID_sorted(i,:)];
                viabilityDataSetForSingleChemical = [viabilityDataSetForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,3:17)];
                concentrationForViabilityForSingleChemical = [concentrationForViabilityForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,18:32)];
            end
            
        else
            averagedResponseData(j,:) = nanmean(responseDataSetForSingleChemical);
            averaged_response_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedResponseData(j,:)];
            
            averagedViabilityData(j,:) = nanmean(viabilityDataSetForSingleChemical);
            averaged_viability_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedViabilityData(j,:)];
            
            averagedConcentrationDataForResponse(j,:) = nanmean(concentrationForResponseForSingleChemical);
            averagedConcentrationDataForResponseWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForResponse(j,:)];
            
            averagedConcentrationDataForViability(j,:) = nanmean(concentrationForViabilityForSingleChemical);
            averagedConcentrationDataForViabilityWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForViability(j,:)];   
            
            j = j+1;
            
            responseDataSetForSingleChemical = [];
            viabilityDataSetForSingleChemical = [];
            concentrationForResponseForSingleChemical = [];
            concentrationForViabilityForSingleChemical = [];
            
            if (contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "antagonist"))
                responseDataForReplicatedExperimentsWithPubIDRankOrder = [responseDataForReplicatedExperimentsWithPubIDRankOrder;...
                    responseDataWithRankAndPubchemSID_sorted(i,:)];
                responseDataSetForSingleChemical = [responseDataSetForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,3:17)];
                concentrationForResponseForSingleChemical = [concentrationForResponseForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,18:32)];
            end
            
            if(contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "viability"))
                viabilityDataForReplicatedExperimentsWithPubIDRankOrder = [viabilityDataForReplicatedExperimentsWithPubIDRankOrder;...
                    responseDataWithRankAndPubchemSID_sorted(i,:)];
                viabilityDataSetForSingleChemical = [viabilityDataSetForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,3:17)];
                concentrationForViabilityForSingleChemical = [concentrationForViabilityForSingleChemical;...
                    responseDataWithRankAndPubchemSID_sorted(i,18:32)];
            end
        end
    else
        if (contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "antagonist"))
            responseDataForReplicatedExperimentsWithPubIDRankOrder = [responseDataForReplicatedExperimentsWithPubIDRankOrder;...
                responseDataWithRankAndPubchemSID_sorted(i,:)];
            responseDataSetForSingleChemical = [responseDataSetForSingleChemical;...
                responseDataWithRankAndPubchemSID_sorted(i,3:17)];
            concentrationForResponseForSingleChemical = [concentrationForResponseForSingleChemical;...
                responseDataWithRankAndPubchemSID_sorted(i,18:32)];
        end
        
        if(contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "viability"))
            viabilityDataForReplicatedExperimentsWithPubIDRankOrder = [viabilityDataForReplicatedExperimentsWithPubIDRankOrder;...
                responseDataWithRankAndPubchemSID_sorted(i,:)];
            viabilityDataSetForSingleChemical = [viabilityDataSetForSingleChemical;...
                responseDataWithRankAndPubchemSID_sorted(i,3:17)];
            concentrationForViabilityForSingleChemical = [concentrationForViabilityForSingleChemical;...
                responseDataWithRankAndPubchemSID_sorted(i,18:32)];
        end
    end
end
 
if (i == length(responseDataWithRankAndPubchemSID_sorted(:,1)) && responseDataWithRankAndPubchemSID_sorted(i,2) == responseDataWithRankAndPubchemSID_sorted(i-1,2))
    averagedResponseData(j,:) = nanmean(responseDataSetForSingleChemical);
    averaged_response_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedResponseData(j,:)];
    
    averagedViabilityData(j,:) = nanmean(viabilityDataSetForSingleChemical);
    averaged_viability_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedViabilityData(j,:)];
    
    averagedConcentrationDataForResponse(j,:) = nanmean(concentrationForResponseForSingleChemical);
    averagedConcentrationDataForResponseWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForResponse(j,:)];
    
    averagedConcentrationDataForViability(j,:) = nanmean(concentrationForViabilityForSingleChemical);
    averagedConcentrationDataForViabilityWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForViability(j,:)];  
else
    
    if (contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "antagonist"))
        averagedResponseData(j,:) = responseDataSetForSingleChemical;
        averaged_response_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedResponseData(j,:)];
        
        averagedConcentrationDataForResponse(j,:) = nanmean(concentrationForResponseForSingleChemical);
        averagedConcentrationDataForResponseWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForResponse(j,:)];
    end
    
    if(contains(sampleDataType(responseDataWithRankAndPubchemSID_sorted(i,1)), "viability"))
        averagedViabilityData(j,:) = viabilityDataSetForSingleChemical;
        averaged_viability_data_with_SID_exclude_concentration(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedViabilityData(j,:)];
        
        averagedConcentrationDataForViability(j,:) = nanmean(concentrationForViabilityForSingleChemical);
        averagedConcentrationDataForViabilityWithSID(j,:) = [responseDataWithRankAndPubchemSID_sorted(i-1,2),averagedConcentrationDataForViability(j,:)];
    end
end
 

fprintf('Row number for Average Response Data Set is : %d.\n', length(averagedResponseData(:,1)))
fprintf('Column number for Average Response Data Set is : %d.\n', length(averagedResponseData(1,:)))
fprintf('Row number for Average Viability Data Set is : %d.\n', length(averagedViabilityData(:,1)))
fprintf('Column number for Average Viability Data Set is : %d.\n', length(averagedViabilityData(1,:)))
 
toc

end