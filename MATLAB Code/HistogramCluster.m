%%
function [] = HistogramCluster(averaged_response_data_with_SID_exclude_concentration,responseDataForReplicatedExperimentsWithPubIDRankOrder,viabilityDataForReplicatedExperimentsWithPubIDRankOrder,averaged_viability_data_with_SID_exclude_concentration)

figure(100)
data_temp = responseDataForReplicatedExperimentsWithPubIDRankOrder(:,2);
data_temp_unique = unique(data_temp);
N = numel(data_temp_unique);
 
 count = zeros(N,1);
 for k = 1:N
     count(k) = sum(data_temp == data_temp_unique(k));
 end
histogram(count)
xlabel("Number of replicates")
ylabel("Frequency")


figure(101)
if (isempty(viabilityDataForReplicatedExperimentsWithPubIDRankOrder)~= 1) % if the data set contains the viability data
data_temp = viabilityDataForReplicatedExperimentsWithPubIDRankOrder(:,2);
data_temp_unique = unique(data_temp);
N = numel(data_temp_unique);
 
 count = zeros(N,1);
 for k = 1:N
     count(k) = sum(data_temp == data_temp_unique(k));
 end
histogram(count)
xlabel("Number of replicates")
ylabel("Frequency")
end

var_vec = [];

for i = 1:1:length(averaged_response_data_with_SID_exclude_concentration(:,1)) 
  obs = averaged_response_data_with_SID_exclude_concentration(i,2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
  var_vec(i) = nanvar(obs); % this is for the variance(n-1), which cannot be too small
end

figure(102)
scatter(randperm(length(var_vec(1,:))),log10(var_vec),'.')
xlabel("Randomized chemical index")
ylabel("Log10 Variance")
hold on

figure(103)
S = quantile(var_vec,[0.1,0.25,0.5,0.75,0.9])
histogram(log10(var_vec),100)
xlabel("Log10 Variance")
ylabel("Frequency")
hold on

figure(104)
data_temp = averaged_response_data_with_SID_exclude_concentration;
data_temp = data_temp(:,2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
data_temp = data_temp(:);
data_temp=data_temp(~isnan(data_temp));
scatter(randperm(length(data_temp)),(data_temp),'.')
xlabel("Randomized chemical-dose index")
ylabel("Averaged response of a chemical at a dose")
hold on

figure(105)
S = quantile(data_temp,[0.1,0.25,0.5,0.75,0.9])
histogram((data_temp),100, "BinWidth", 1)
xlabel("Averaged response of a chemical at a dose")
ylabel("Count")
hold on

Max_res = [];
for i = 1:1:length(averaged_response_data_with_SID_exclude_concentration(:,1))
    averaged_response_data_temp = averaged_response_data_with_SID_exclude_concentration(:,2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
    Max_res = [Max_res,max(averaged_response_data_temp(i,:))];
end

figure(106)
scatter(Max_res,log10(var_vec),'.')
xlabel("Maximal response")
ylabel("Log10 Variance")
% xline(15,'r')
% yline(1.65,'r')
hold on

Min_res = [];
for i = 1:1:length(averaged_response_data_with_SID_exclude_concentration(:,1))
    averaged_response_data_temp = averaged_response_data_with_SID_exclude_concentration(:,2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
    Min_res = [Min_res,min(averaged_response_data_temp(i,:))];
end

figure(107)
scatter(Min_res,log10(var_vec),'.')
xlabel("Minimal response")
ylabel("Log10 Variance")
% xline(15,'r')
% yline(1.65,'r')
hold on

figure(108)
scatter(Min_res,Max_res,'.')
xlabel("Minimal response")
ylabel("Maximal response")
% xline(15,'r')
% yline(1.65,'r')
hold on

figure(109)
data_temp = averaged_viability_data_with_SID_exclude_concentration;
data_temp = data_temp(:,2:length(averaged_viability_data_with_SID_exclude_concentration(1,:)));
data_temp = data_temp(:);
data_temp=data_temp(~isnan(data_temp));
scatter(randperm(length(data_temp)),(data_temp),'.')
xlabel("Randomized chemical-dose index")
ylabel("Averaged viability of a chemical at a dose")
hold on

figure(120)
S = quantile(data_temp,[0.1,0.25,0.5,0.75,0.9])
histogram((data_temp),100, "BinWidth", 1)
xlabel("Averaged viability of a chemical at a dose")
ylabel("Count")
hold on

end