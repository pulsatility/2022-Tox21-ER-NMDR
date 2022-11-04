% MatLab code for "Identification of nonmonotonic concentration-responses
% in Tox21 high- throughput screening estrogen receptor assays" published
% in Toxicology and Applied Pharmcology 2022

clear all
clc
close all

%% ---Select assay number to start--- %%
assayNumber = 8;


%% Threshold values used for filtering and classification

number_of_concens_response = 15;
U_shape_magnitude = 15;
Bell_shape_magnitude = 15;
Flat_magnitude = 30;
interference_boundary = -10;
autofluorescence_boundary = 10;
viability_correlation_threshold = 0.4;
interference_correlation_threshold = 0.4;
interference_vertical_drop_threshold = 1.5*abs(interference_boundary);


%% Import all Tox21 ER and related data and average activities of replicates 
%Read the Tripod txt file. The files can be downloaded from https://tripod.nih.gov/tox.

% assay #1
Data_strings_er_luc_agonist_p4 = fileread('tox21-er-luc-bg1-4e2-agonist-p4.txt');
Data_cells_er_luc_agonist_p4 = textscan(Data_strings_er_luc_agonist_p4, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_luc_agonist_p4,full_response_data_er_luc_agonist_p4,...
    averaged_viability_data_er_luc_agonist_p4,full_viability_data_er_luc_agonist_p4,...
    averaged_concentration_response_er_luc_agonist_p4,...
    averaged_concentration_viability_er_luc_agonist_p4] = CalculateAverageDataResponse3(Data_cells_er_luc_agonist_p4);

% assay #2 
Data_strings_er_luc_agonist_p2 = fileread('tox21-er-luc-bg1-4e2-agonist-p2.txt');
Data_cells_er_luc_agonist_p2 = textscan(Data_strings_er_luc_agonist_p2, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_luc_agonist_p2,full_response_data_er_luc_agonist_p2,...
    averaged_viability_data_er_luc_agonist_p2,full_viability_data_er_luc_agonist_p2,...
    averaged_concentration_response_er_luc_agonist_p2,...
    averaged_concentration_viability_er_luc_agonist_p2] = CalculateAverageDataResponse2(Data_cells_er_luc_agonist_p2);

% assay #3
Data_strings_er_luc_antagonist_p2 = fileread('tox21-er-luc-bg1-4e2-antagonist-p2.txt');
Data_cells_er_luc_antagonist_p2 = textscan(Data_strings_er_luc_antagonist_p2, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_luc_antagonist_p2,full_response_data_er_luc_antagonist_p2,...
    averaged_viability_data_er_luc_antagonist_p2,full_viability_data_er_luc_antagonist_p2,...
    averaged_concentration_response_er_luc_antagonist_p2,...
    averaged_concentration_viability_er_luc_antagonist_p2]= CalculateAverageDataResponse4(Data_cells_er_luc_antagonist_p2);

% assay #4
Data_strings_er_luc_antagonist_p1 = fileread('tox21-er-luc-bg1-4e2-antagonist-p1.txt');
Data_cells_er_luc_antagonist_p1 = textscan(Data_strings_er_luc_antagonist_p1, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_luc_antagonist_p1,full_response_data_er_luc_antagonist_p1,...
    averaged_viability_data_er_luc_antagonist_p1,full_viability_data_er_luc_antagonist_p1,...
    averaged_concentration_response_er_luc_antagonist_p1,...
    averaged_concentration_viability_er_luc_antagonist_p1]= CalculateAverageDataResponse4(Data_cells_er_luc_antagonist_p1);

% assay #5
Data_strings_er_bla_agonist_p2 = fileread('tox21-er-bla-agonist-p2.txt');
Data_cells_er_bla_agonist_p2 = textscan(Data_strings_er_bla_agonist_p2, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_bla_agonist_p2,full_response_data_er_bla_agonist_p2,...
    averaged_viability_data_er_bla_agonist_p2,full_viability_data_er_bla_agonist_p2,...
    averaged_concentration_response_er_bla_agonist_p2,...
    averaged_concentration_viability_er_bla_agonist_p2] = CalculateAverageDataResponse1(Data_cells_er_bla_agonist_p2);

% assay #6
Data_strings_er_bla_antagonist_p1 = fileread('tox21-er-bla-antagonist-p1.txt');
Data_cells_er_bla_antagonist_p1 = textscan(Data_strings_er_bla_antagonist_p1, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_er_bla_antagonist_p1,full_response_data_er_bla_antagonist_p1,...
    averaged_viability_data_er_bla_antagonist_p1,full_viability_data_er_bla_antagonist_p1,...
    averaged_concentration_response_er_bla_antagonist_p1,...
    averaged_concentration_viability_er_bla_antagonist_p1] = CalculateAverageDataResponse1(Data_cells_er_bla_antagonist_p1);

% assay #7
Data_strings_erb_bla_p1 = fileread('tox21-erb-bla-p1.txt');
Data_cells_erb_bla_p1 = textscan(Data_strings_erb_bla_p1, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_erb_bla_p1,full_response_data_erb_bla_p1,...
    averaged_viability_data_erb_bla_p1,full_viability_data_erb_bla_p1,...
    averaged_concentration_response_erb_bla_p1,...
    averaged_concentration_viability_erb_bla_p1] = CalculateAverageDataResponse1(Data_cells_erb_bla_p1);

% assay #8
Data_strings_erb_bla_antagonist_p1 = fileread('tox21-erb-bla-antagonist-p1.txt');
%Convert the strings to MatLab Cells: %s is for string (characters) data type, %f is for numerical (double) data type
Data_cells_erb_bla_antagonist_p1 = textscan(Data_strings_erb_bla_antagonist_p1, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_response_data_erb_bla_antagonist_p1,full_response_data_erb_bla_antagonist_p1,...
    averaged_viability_data_erb_bla_antagonist_p1,full_viability_data_erb_bla_antagonist_p1,...
    averaged_concentration_response_erb_bla_antagonist_p1,...
    averaged_concentration_viability_erb_bla_antagonist_p1] = CalculateAverageDataResponse1(Data_cells_erb_bla_antagonist_p1);

% interference data
Tox21_interference_data = fileread('tox21-luc-biochem-p1.txt');
Data_cells_Tox21_interference_data = textscan(Tox21_interference_data, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_interference_data, averaged_concentration_for_interference, full_interference_data] = CalculateAverageDataResponse5(Data_cells_Tox21_interference_data);

% autofluorescence
Tox21_autofluorescence_data = fileread('tox21-spec-hek293-p1.txt');
Data_cells_Tox21_autofluorescence_data = textscan(Tox21_autofluorescence_data, '%s%s%s%s%s%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%s%s%s%s%s%s','Delimiter', '\t', 'Headerlines', 1 );
[averaged_autofluorescence_data, averaged_concentration_for_autofluorescence, full_autofluorescence_data] = CalculateAverageDataResponse6(Data_cells_Tox21_autofluorescence_data);


%% Convert specific assay names to common variables

switch assayNumber
    
    case assayNumber == 1
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -10;
        viability_upper_bound = 10;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_agonist_p4;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_agonist_p4;
        
        averaged_concentration_response = averaged_concentration_response_er_luc_agonist_p4;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_agonist_p4;
        
        full_response_data = full_response_data_er_luc_agonist_p4;
        full_viability_data = full_viability_data_er_luc_agonist_p4;
    
        
    case assayNumber == 2
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -10;
        viability_upper_bound = 10;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_agonist_p2;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p1; % use assay 4 viability data
        
        averaged_concentration_response = averaged_concentration_response_er_luc_agonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p1;
        
        full_response_data = full_response_data_er_luc_agonist_p2;
        full_viability_data = full_viability_data_er_luc_antagonist_p1;
        
        
    case assayNumber == 3
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -10;
        viability_upper_bound = 10;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_antagonist_p2;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p2;
        
        averaged_concentration_response = averaged_concentration_response_er_luc_antagonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p2;
        
        full_response_data = full_response_data_er_luc_antagonist_p2;
        full_viability_data = full_viability_data_er_luc_antagonist_p2;
    
        
    case assayNumber == 4
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -10;
        viability_upper_bound = 10;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p1;
        
        averaged_concentration_response = averaged_concentration_response_er_luc_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p1;
        
        full_response_data = full_response_data_er_luc_antagonist_p1;
        full_viability_data = full_viability_data_er_luc_antagonist_p1;
        
        
    case assayNumber == 5
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -20;
        viability_upper_bound = 20;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_bla_agonist_p2;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_bla_antagonist_p1; % use assay 6 viability data
        
        averaged_concentration_response = averaged_concentration_response_er_bla_agonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_bla_antagonist_p1;
        
        full_response_data = full_response_data_er_bla_agonist_p2;
        full_viability_data = full_viability_data_er_bla_antagonist_p1;

        
    case assayNumber == 6
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -20;
        viability_upper_bound = 20;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_bla_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_bla_antagonist_p1;
        
        averaged_concentration_response = averaged_concentration_response_er_bla_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_er_bla_antagonist_p1;
        
        full_response_data = full_response_data_er_bla_antagonist_p1;
        full_viability_data = full_viability_data_er_bla_antagonist_p1;
        
        
    case assayNumber == 7
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -15;
        viability_upper_bound = 15;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_erb_bla_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_erb_bla_p1;
        
        averaged_concentration_response = averaged_concentration_response_erb_bla_p1;
        averaged_concentration_viability = averaged_concentration_viability_erb_bla_p1;
        
        full_response_data = full_response_data_erb_bla_p1;
        full_viability_data = full_viability_data_erb_bla_p1;

        
    case assayNumber == 8
        
        % cutoff values
        upper_boundary = 20;
        lower_boundary = -20;
        magnitude_of_outlier = 1.5*upper_boundary;
        
        viability_lower_bound = -15;
        viability_upper_bound = 15;
        viability_vertical_drop_threshold = 1.5*viability_upper_bound;
        
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        
        % data input
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_erb_bla_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_erb_bla_antagonist_p1;
        
        averaged_concentration_response = averaged_concentration_response_erb_bla_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_erb_bla_antagonist_p1;
        
        full_response_data = full_response_data_erb_bla_antagonist_p1;
        full_viability_data = full_viability_data_erb_bla_antagonist_p1;

end


%% Plot histograms and scatter plot data points

HistogramCluster(averaged_response_data_with_SID_exclude_concentration,full_response_data,full_viability_data,averaged_viability_data_with_SID_exclude_concentration);


figure(1000)
data_temp = averaged_interference_data;
data_temp = data_temp(:,2:length(averaged_interference_data(1,:)));
data_temp = data_temp(:);
data_temp=data_temp(~isnan(data_temp));
scatter(randperm(length(data_temp)),(data_temp),'.')
xlabel("Randomized chemical-dose index")
ylabel("Averaged interference of a chemical at a dose")
hold on

figure(1010)
S = quantile(data_temp,[0.1,0.25,0.5,0.75,0.9])
histogram((data_temp),100, "BinWidth", 1)
xlabel("Averaged interference of a chemical at a dose")
ylabel("Count")
hold on

figure(2000)
data_temp = averaged_autofluorescence_data;
data_temp = data_temp(:,2:length(averaged_autofluorescence_data(1,:)));
data_temp = data_temp(:);
data_temp=data_temp(~isnan(data_temp));
scatter(randperm(length(data_temp)),(data_temp),'.')
xlabel("Randomized chemical-dose index")
ylabel("Averaged autofluorescence of a chemical at a dose")
hold on

figure(2010)
S = quantile(data_temp,[0.1,0.25,0.5,0.75,0.9])
histogram((data_temp),100, "BinWidth", 1)
xlabel("Averaged autofluorescence of a chemical at a dose")
ylabel("Count")
hold on


%% (A) Initial screening 

var_vec = [];
selected_LowVariance_pubchemSID = [];
for i = 1:1:length(averaged_response_data_with_SID_exclude_concentration(:,1))
  obs = averaged_response_data_with_SID_exclude_concentration(i,2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
  var_vec(i) = nanvar(obs);
    
  if (max(obs) <= upper_boundary && min(obs) >= lower_boundary)
  
    selected_LowVariance_pubchemSID = [selected_LowVariance_pubchemSID, averaged_response_data_with_SID_exclude_concentration(i,1)];
  end
end

fprintf('Number of Screened Out Low-Variance Response Curves is : %d.\n', length(selected_LowVariance_pubchemSID));
m = length(averaged_response_data_with_SID_exclude_concentration(:,1)) - length(selected_LowVariance_pubchemSID);
fprintf('Number of Remaining Response Curves is : %d.\n', m);

averaged_response_data_exclude_low_variance = averaged_response_data_with_SID_exclude_concentration;
indset = [];
for i = 1:1:length(selected_LowVariance_pubchemSID)
    indset(i) = find(averaged_response_data_with_SID_exclude_concentration(:,1) == selected_LowVariance_pubchemSID(i));
end

averaged_response_data_exclude_low_variance(indset,:) = [];

% % View the concentration-response curves for screened out low-activity chemicals
% for i = 1:1:length(selected_LowVariance_pubchemSID)
%     seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == selected_LowVariance_pubchemSID(i));
%     single_response_data = averaged_response_data_with_SID_exclude_concentration(seq_in_response_data,...
%         2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
%     conc = averaged_concentration_response(seq_in_response_data,...
%         2:length(averaged_concentration_response(1,:)));
%     plot(log10(conc),single_response_data,'-ko','LineWidth',2)
%     %xlim([-4 3])
%     ylim([-100 100])
%     %yticks([-100:20:100])
%     xlabel('Dose')
%     ylabel('Response')
%     set(gca,'fontsize',18)
%     box off
%     grid on
%     hold on
% end


% Screening out response curves with outlier spike

selected_HighVariance_pubchemSID = [];

for i = 1:1:length(averaged_response_data_exclude_low_variance(:,1))
    flag = 0;
    obs = averaged_response_data_exclude_low_variance(i,2:length(averaged_response_data_exclude_low_variance(1,:)));
    obs = obs(~isnan(obs));
    obs_diff = diff(obs);
    
    maximal_point_index = find(obs == max(obs));
    obs_exclude_maximal = obs;
    obs_exclude_maximal(maximal_point_index)=[];
    second_maximal_point_index = find(obs_exclude_maximal == max(obs_exclude_maximal));
    
    minimal_point_index = find(obs == min(obs));
    obs_exclude_minimal = obs;
    obs_exclude_minimal(minimal_point_index)=[];
    second_minimal_point_index = find(obs_exclude_minimal == min(obs_exclude_minimal));
    
    
    for j = 2:1:length(obs_diff)
        % the sign in the change between the two neighor points has to be different
        % if the change happens at the last point; don't count
        % don't count twice for one observation; this is the reason why a "flag" variable is set
        if ((obs_diff(j)*obs_diff(j-1) < 0 && ...
                abs(obs_diff(j-1))>= magnitude_of_outlier &&...
                abs(obs_diff(j))>= magnitude_of_outlier &&...
                maximal_point_index == j &&...% the outlier needs to be the maximal or minimal point
                abs(obs(1)-obs(end)) <= 100 &&... % exclude the elevated curves
                abs(obs_exclude_maximal(second_maximal_point_index)-obs_exclude_maximal(1))...
                <= 100 &&... % exclude bell shape curves
                j ~= length(obs_diff)) ||...
                (obs_diff(j)*obs_diff(j-1) < 0 && ...
                abs(obs_diff(j-1))>= magnitude_of_outlier &&...
                abs(obs_diff(j))>= magnitude_of_outlier &&...
                minimal_point_index == j &&...% the outlier needs to be the maximal or minimal point
                abs(obs(1)-obs(end)) <= 100 &&... % exclude the elevated curves
                abs(obs_exclude_minimal(second_minimal_point_index)-obs_exclude_minimal(1))...
                <= 100 &&... % exclude U shape curves
                j ~= length(obs_diff)))
            if (flag == 0)
                selected_HighVariance_pubchemSID = [selected_HighVariance_pubchemSID,averaged_response_data_exclude_low_variance(i,1)];
                flag = 1;
            end
        end
    end
end

fprintf('Number of Screened Out Response Curves with One Outstanding Peak is : %d.\n', length(selected_HighVariance_pubchemSID));
m = length(averaged_response_data_exclude_low_variance(:,1)) - length(selected_HighVariance_pubchemSID);
fprintf('Remaining Response Curves after Screening Out Low-variance and Outstanding Peak is : %d.\n', m);

averaged_response_data_exclude_low_variance_outstanding_peak = averaged_response_data_exclude_low_variance;
indset = [];
for i = 1:1:length(selected_HighVariance_pubchemSID)
    indset(i) = find(averaged_response_data_exclude_low_variance(:,1) == selected_HighVariance_pubchemSID(i));
end

averaged_response_data_exclude_low_variance_outstanding_peak(indset,:) = [];

% View the concentration-response curves for screened out single-data point peak chemicals
for i = 1:1:length(selected_HighVariance_pubchemSID)
   figure(i)
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        selected_HighVariance_pubchemSID(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        selected_HighVariance_pubchemSID(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            selected_HighVariance_pubchemSID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            selected_HighVariance_pubchemSID(i));
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == selected_HighVariance_pubchemSID(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == selected_HighVariance_pubchemSID(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == selected_HighVariance_pubchemSID(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == selected_HighVariance_pubchemSID(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data/sqrt(length(chemical_index_in_full_response_data)),'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data/sqrt(length(chemical_index_in_full_viability_data)),'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data/sqrt(length(chemical_index_in_full_interference_data)),'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data/sqrt(length(chemical_index_in_full_autofluoresence_data)),'-bo','LineWidth',2)
    end
    
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(selected_HighVariance_pubchemSID(i))],'Location','southwest');
    legend('boxoff')
end
    

%% Obtain moving average to reduce the local effect

mov_para_1 = 5;
mov_para_2 = 2;

AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID =...
averaged_response_data_exclude_low_variance_outstanding_peak(:,2:length(averaged_response_data_exclude_low_variance_outstanding_peak(1,:)));
AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1 = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID;
AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2 = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID;

for i = 1:1:length(AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(:,1))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(i,:);
    single_response_exclude_nan = single_response_data;
    single_response_exclude_nan=single_response_exclude_nan(~isnan(single_response_exclude_nan));
    single_response_MA_1 = movmean(single_response_exclude_nan,mov_para_1);
    single_response_MA_2 = movmean(single_response_exclude_nan,mov_para_2);
    k = 1;
    for j = 1:1:number_of_concens_response
        if (~isnan(single_response_data(j)))
        AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(i,j) = single_response_MA_1(k);
        AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(i,j) = single_response_MA_2(k);
        k = k + 1;    
        end
    end    
end

% View the comparison between orignial plots and plots after applying the moving Average
start_of_obs = 1;
end_of_obs = length(AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(:,1));

for i = 171:1:172
    figure(i)
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        averaged_response_data_exclude_low_variance_outstanding_peak(i,1));
    single_response_data = averaged_response_data_with_SID_exclude_concentration(seq_in_response_data,...
        2:length(averaged_response_data_with_SID_exclude_concentration(1,:)));
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));

    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(i,:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(i,:);
    plot(log10(conc), single_response_data,'-ko','LineWidth',2)
    hold on
    plot(log10(conc),single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(log10(conc),single_response_MA_2,'-go','LineWidth',2)
    legend('orignial',['movavg =' num2str(mov_para_1)],['movavg =' num2str(mov_para_2)],...
        'Location','northeast')
    %xlim([-4 3])
    %ylim([-100 100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on   
end


%% (B) Recursive, correlation-based clustering

Cluster_ID_Matrix = [];
current_chemical_set = [1:1:length(AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(:,1))];
k = 1;

tic

while (length(current_chemical_set) > 0)
    selected_index_temp_2 = 1;
    % initialize corr_temp
    corr_temp = 0;
    % Step 1: Pick a random chemical    
    % pick a random chemical
    % randsample(n,k) returns k values sampled uniformly at random, without replacement, from the integers 1 to n.
    % note: a sampled index is not equal to an actual index.
    % For example, the actual chemical index in a chemical set [1,2,3,5] is 5 when the 4th one is picked.
    sampled_index = randsample(length(current_chemical_set), 1);
    actual_index = current_chemical_set(sampled_index);
    current_chemical_set_exclude_chosen_one = current_chemical_set;
    current_chemical_set_exclude_chosen_one(sampled_index)=[];
    Cluster_ID_Matrix(k,1) = actual_index;
    position_in_column = 2;

    % Step 2: Loop through the chemical set that exclude the chose one, and
    % pikc the one that has the largest correlation with it

    for i = 1:1:length(current_chemical_set_exclude_chosen_one)
        current_dose_response_curve = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(actual_index,:);
        compared_dose_response_curve = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(current_chemical_set_exclude_chosen_one(i),:);
        
        current_dose_response_curve = current_dose_response_curve(~isnan(current_dose_response_curve));
        compared_dose_response_curve = compared_dose_response_curve(~isnan(compared_dose_response_curve));
        
        % the length of current dose response curve may be not as the same
        % as the compared ones; so compare the response at higher
        % concentrations
        if (length(current_dose_response_curve)>length(compared_dose_response_curve))
            current_dose_response_curve = current_dose_response_curve((end-length(compared_dose_response_curve)+1):end);
        else
            compared_dose_response_curve = compared_dose_response_curve((end-length(current_dose_response_curve)+1):end);
        end
        
        if (corr(current_dose_response_curve',compared_dose_response_curve') > corr_temp)
            corr_temp = corr(current_dose_response_curve',compared_dose_response_curve');
            selected_index_temp_1 = current_chemical_set_exclude_chosen_one(i);
        end  
    end
    
    if (length(current_chemical_set_exclude_chosen_one) > 0 && corr_temp > 0.75)
        % Store the selected chemical into a cluster set
        Cluster_ID_Matrix(k,position_in_column) = selected_index_temp_1;
        position_in_column = position_in_column + 1;
    end
       
    % Step 3: caclulate the correlation between each chemical in the existing 
    % clustering set and chemicals from chemical set that excludes the existing cluster
    % For each chemical in the chemical set that excludes the existing
    % cluster, compute correlation between the chemical and the chemicals in 
    % the existing chemical set, find a chemimcal whose minimal correlation is maximal 
    % among all the chemicals in the chemical set that excludes the exisiting cluster
    % repeat this process until the maximal-minial correlation is less than
    % a predefined threshold

    if (length(current_chemical_set_exclude_chosen_one) > 1)
        % initialize corr_temp_2 that use to store the maximal correlation
        corr_temp_2 = 1;
        while (corr_temp_2 > 0.75)
            Cluster_ID_Matrix_vector = Cluster_ID_Matrix(k,:);
            Cluster_ID_Matrix_vector(Cluster_ID_Matrix_vector == 0) = [];
            Chemical_Set_Exclude_Selected_Ones = current_chemical_set;
            % Exclude the selected chemicals from exisiting chemical set
            Cluster_ID_Matrix_vector_Index = [];
            for i = 1:1:length(Cluster_ID_Matrix_vector)
                Cluster_ID_Matrix_vector_Index = find(Chemical_Set_Exclude_Selected_Ones == Cluster_ID_Matrix_vector(i));
                Chemical_Set_Exclude_Selected_Ones (Cluster_ID_Matrix_vector_Index) = [];
            end
            
            corr_temp_1 = 1;
            corr_temp_2 = 0;
            for j = 1:1:length(Chemical_Set_Exclude_Selected_Ones)
                % initialize corr_temp_1 for each chemical
                corr_temp_1 = 1;
                for m = 1:1:length(Cluster_ID_Matrix_vector)
                    chemical_from_excluded_set = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Chemical_Set_Exclude_Selected_Ones(j),:);
                    chemical_from_excluded_set = chemical_from_excluded_set(~isnan(chemical_from_excluded_set));
                    chemcial_from_exisiting_cluster = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Cluster_ID_Matrix(k,m),:);
                    chemcial_from_exisiting_cluster = chemcial_from_exisiting_cluster(~isnan(chemcial_from_exisiting_cluster));
                    
                    % the length of current dose response curve may be not as the same
                    % as the compared ones; so compare the response at higher
                    % concentrations
                    if (length(chemical_from_excluded_set)>length(chemcial_from_exisiting_cluster))
                        chemical_from_excluded_set = chemical_from_excluded_set((end-length(chemcial_from_exisiting_cluster)+1):end);
                    else
                        chemcial_from_exisiting_cluster = chemcial_from_exisiting_cluster((end-length(chemical_from_excluded_set)+1):end);
                    end
                    
                    % corr_temp_1 store the minimal correlation
                    if (corr(chemical_from_excluded_set',chemcial_from_exisiting_cluster') < corr_temp_1)
                        corr_temp_1 = corr(chemical_from_excluded_set',chemcial_from_exisiting_cluster');
                    end
                end
                
                % corr_temp_2 store the maixmal correlation
                if (corr_temp_1 > corr_temp_2)
                    corr_temp_2 = corr_temp_1;
                    selected_index_temp_2 = Chemical_Set_Exclude_Selected_Ones(j);
                end
            end
            
            if(corr_temp_2 > 0.75)
                Cluster_ID_Matrix(k,position_in_column) = selected_index_temp_2;
                position_in_column = position_in_column + 1;
            end
            
            corr_temp_2
        end
    end
    
    % Step 4: Exclude the selected chemicals from the chemical set and
    % Repeat Step 1
    Cluster_ID_Matrix_vector =Cluster_ID_Matrix(k,:);
    Cluster_ID_Matrix_vector(Cluster_ID_Matrix_vector == 0) = [];
    Cluster_ID_Matrix_vector_Index = [];
    for i = 1:1:length(Cluster_ID_Matrix_vector)
        Cluster_ID_Matrix_vector_Index = find(current_chemical_set == Cluster_ID_Matrix_vector(i));
        current_chemical_set (Cluster_ID_Matrix_vector_Index) = [];
    end
    
    k = k + 1;  
end 

toc
  
number_of_clusters = Cluster_ID_Matrix(:,1);
number_of_clusters = length(number_of_clusters);
fprintf('Number of Clusters is : %d.\n', number_of_clusters);


% Plot CRCs based on clustering

for i = 1:1:number_of_clusters
    figure(i)
    Cluster_ID_Set = Cluster_ID_Matrix(i,:);
    Cluster_ID_Set(Cluster_ID_Set == 0) = [];
    
    for j = 1:1:length(Cluster_ID_Set)
        response_for_one_chemical = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID...
            (Cluster_ID_Set(j),:);
        
        seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
            averaged_response_data_exclude_low_variance_outstanding_peak(Cluster_ID_Set(j),1));
        
        conc = averaged_concentration_response(seq_in_response_data,...
            2:length(averaged_concentration_response(1,:)));
        
        plot(log10(conc),response_for_one_chemical,'-ko','LineWidth',2)
        %xlim([-4 3])
%         ylim([-100 100])
%         yticks([-100:20:100])
        xlabel('Dose')
        ylabel('Response')
        set(gca,'fontsize',18)
        box off
        grid on
        hold on
    end
    
end

% After visual insepction, the following six shapes for CRCs were identified:
% Flat, monotonic decrease, monotonic increase, Bell, U, and S curves


%% (C) Pattern-restricted classification of CRCs 

Flat_index = [];
Increasing_index = [];
Decreasing_index = [];
U_index_mov1 = [];
U_index_mov2 = [];
U_index = [];
Bell_index_mov1 = [];
Bell_index_mov2 = [];
Bell_index = [];
S_index = [];

start_index = 1;
end_index = length(averaged_response_data_exclude_low_variance_outstanding_peak(:,1));


for i = start_index:1:end_index

    flag_1 = 0;
    flag_2 = 0;
    flag_3 = 0;
    flag_4 = 0;
    
    single_response = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(i,:);
    single_response_exclude_NA = single_response;
    single_response_exclude_NA = single_response_exclude_NA(~isnan(single_response_exclude_NA));
    min_index_data = find(single_response_exclude_NA == min(single_response_exclude_NA));
    max_index_data = find(single_response_exclude_NA == max(single_response_exclude_NA));
    
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(i,:);
    single_response_MA_1_exclude_NA = single_response_MA_1;
    single_response_MA_1_exclude_NA = single_response_MA_1_exclude_NA(~isnan(single_response_MA_1_exclude_NA));
    
    first_index = 1;
    last_index = length(single_response_MA_1_exclude_NA);
    
    min_index = find(single_response_MA_1_exclude_NA == min(single_response_MA_1_exclude_NA));
    max_index = find(single_response_MA_1_exclude_NA == max(single_response_MA_1_exclude_NA));
    
   
    % Bell Shape
    if (max_index_data > 2 && max_index_data < last_index - 1 &&...
            single_response_MA_1_exclude_NA(max_index) - single_response_MA_1_exclude_NA(first_index) >...
            Bell_shape_magnitude && single_response_MA_1_exclude_NA(max_index) -...
            single_response_MA_1_exclude_NA(last_index) > Bell_shape_magnitude)
        
            Bell_index_mov1 = [Bell_index_mov1,i];  
            flag_1 = 1;
    end 
     
    % U Shape
    if (min_index_data > 2 && min_index_data < last_index - 1 &&...
          abs(single_response_MA_1_exclude_NA(min_index) - single_response_MA_1_exclude_NA(first_index))...
          > U_shape_magnitude && abs(single_response_MA_1_exclude_NA(min_index) -...
          single_response_MA_1_exclude_NA(last_index)) > U_shape_magnitude)

            U_index_mov1 = [U_index_mov1,i];  
            flag_2 = 1;
    end
    
    if (flag_1 == 0 && flag_2 == 0)

        single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(i,:);
        single_response_MA_2_exclude_NA = single_response_MA_2;
        single_response_MA_2_exclude_NA = single_response_MA_2_exclude_NA(~isnan(single_response_MA_2_exclude_NA));
        
        first_index = 1;
        last_index = length(single_response_MA_2_exclude_NA);
        
        min_index = find(single_response_MA_2_exclude_NA == min(single_response_MA_2_exclude_NA));
        max_index = find(single_response_MA_2_exclude_NA == max(single_response_MA_2_exclude_NA));
        
      
        % Bell Shape
        if (max_index_data > 2 && max_index_data < last_index - 1 &&...
                single_response_MA_2_exclude_NA(max_index) - single_response_MA_2_exclude_NA(first_index) >...
                Bell_shape_magnitude && single_response_MA_2_exclude_NA(max_index) -...
                single_response_MA_2_exclude_NA(last_index) > Bell_shape_magnitude)
            
            Bell_index_mov2 = [Bell_index_mov2,i];
            flag_3 = 1;
        end
        
             
        % U Shape
        if (min_index_data > 2 && min_index_data < last_index - 1 &&...
                abs(single_response_MA_2_exclude_NA(min_index) - single_response_MA_2_exclude_NA(first_index))...
                > U_shape_magnitude && abs(single_response_MA_2_exclude_NA(min_index) -...
                single_response_MA_2_exclude_NA(last_index)) > U_shape_magnitude)
            
            U_index_mov2 = [U_index_mov2,i];
            flag_4 = 1;
        end
        
        if (flag_3 == 0 && flag_4 == 0)
            if (abs(single_response_exclude_NA(max_index_data) - single_response_exclude_NA(min_index_data))<= Flat_magnitude)
                Flat_index = [Flat_index,i];
            else
                if (min_index_data < max_index_data)
                    Increasing_index = [Increasing_index,i];
                else
                    Decreasing_index = [Decreasing_index,i];
                end
            end  
        end
        
        
        if (flag_3 == 1 && flag_4 == 1)   
            S_index = [S_index,i];
%             U_index(find(U_index == i)) = []; % Exclude S curves from U shape index
%             Bell_index(find(Bell_index == i)) = []; % Exclude S curves from Bell shape index       
        end  
    end

    % need to revisit it later
    if (flag_1 == 1 && flag_2 == 1)
       S_index = [S_index,i]; 
%        U_index(find(U_index == i)) = []; % Exclude S curves from U shape index
%        Bell_index(find(Bell_index == i)) = []; % Exclude S curves from Bell shape index     
    end
    
end

% Exclude S curves from U shape index
a=setdiff(U_index_mov1,S_index);
index=find(ismember(U_index_mov1,a));
U_index_mov1=U_index_mov1(index);

a=setdiff(U_index_mov2,S_index);
index=find(ismember(U_index_mov2,a));
U_index_mov2=U_index_mov2(index);

% Exclude S curves from Bell shape index   
a=setdiff(Bell_index_mov1,S_index);
index=find(ismember(Bell_index_mov1,a));
Bell_index_mov1=Bell_index_mov1(index);

a=setdiff(Bell_index_mov2,S_index);
index=find(ismember(Bell_index_mov2,a));
Bell_index_mov2=Bell_index_mov2(index);

U_index = [U_index_mov1, U_index_mov2];
Bell_index = [Bell_index_mov1, Bell_index_mov2];

fprintf('The number of flat shape is : %d.\n', length(Flat_index));
fprintf('The number of decreasing shape is : %d.\n', length(Decreasing_index));
fprintf('The number of increasing shape is : %d.\n', length(Increasing_index));
fprintf('The number of U shape is : %d.\n', length(U_index));
fprintf('The number of Bell shape is : %d.\n', length(Bell_index));
fprintf('The number of S shape is : %d.\n', length(S_index));
fprintf('The number of shape is : %d.\n', length(Flat_index)+length(Decreasing_index)+...
    length(Increasing_index)+length(U_index)+length(Bell_index)+length(S_index));


%% Plot CRCs of different shapes

flat_pubchem_ID_set = [];
increasing_pubchem_ID_set = [];
decreasing_pubchem_ID_set = [];
U_pubchem_ID_set = [];
Bell_pubchem_ID_set = [];
S_pubchem_ID_set = [];

% Flat shape
for i = 1:1:length(Flat_index)
 
    figure(Flat_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(Flat_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Flat_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(Flat_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Flat_index(i),1);
    flat_pubchem_ID_set = [flat_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end

% Increasing shape
for i = 1:1:50%length(Increasing_index)
 
    figure(Increasing_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(Increasing_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Increasing_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(Increasing_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Increasing_index(i),1);
    increasing_pubchem_ID_set = [increasing_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end

% Decrease shape
for i = 1:1:50%length(Decreasing_index)
 
    figure(Decreasing_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(Decreasing_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Decreasing_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(Decreasing_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Decreasing_index(i),1);
    decreasing_pubchem_ID_set = [decreasing_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end

% U shape

for i = 1:1:length(U_index)
 
    figure(U_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(U_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(U_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(U_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(U_index(i),1);
    U_pubchem_ID_set = [U_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    
%     plot(conc,single_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end

% Bell shape
for i = 1:1:50%length(Bell_index)
 
    figure(Bell_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(Bell_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(Bell_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(Bell_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Bell_index(i),1);
    Bell_pubchem_ID_set = [Bell_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
%     plot(conc,single_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end

% S shape
for i = 1:1:length(S_index)
 
    figure(S_index(i))
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(S_index(i),:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(S_index(i),:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(S_index(i),:);
    
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(S_index(i),1);
    S_pubchem_ID_set = [S_pubchem_ID_set,pubchem_ID];
    % find the chemical index in the original averaged response data; the index in
    % the post-processed averaged response data excludes some chemicals
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == pubchem_ID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);

    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    seq_in_response_data = find(averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        pubchem_ID);
    conc = averaged_concentration_response(seq_in_response_data,...
        2:length(averaged_concentration_response(1,:)));
    conc = log10(conc);
        
    errorbar(conc,single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    plot(conc,single_response_MA_1,'-bo','LineWidth',2)
    hold on
    plot(conc,single_response_MA_2,'-go','LineWidth',2)
    hold on
    
%     xlim([-4 3])
%     ylim([-100 120])
%     yticks([-100:20:120])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on  
    legend(['ID=' num2str(pubchem_ID)],'Location','southwest');
    legend('boxoff')
end


%% (D) False-positive NMDR exclusion:(1) Identify compounds active in cytotoxicity, luciferase inhibition, autofluoresecence assays

% Identify chemicals active in cytotoxicity assay (cutoff specified based on assay)
if (length(averaged_viability_data_with_SID_exclude_concentration(1,:)) > 2)
    selected_LowViability_pubchemSID = [];
    for i = 1:1:length(averaged_viability_data_with_SID_exclude_concentration(:,1))
        obs = averaged_viability_data_with_SID_exclude_concentration(i,2:length(averaged_viability_data_with_SID_exclude_concentration(1,:)));
        
        if (max(obs) <= viability_upper_bound && min(obs) >= viability_lower_bound)
            
            selected_LowViability_pubchemSID = [selected_LowViability_pubchemSID, averaged_viability_data_with_SID_exclude_concentration(i,1)];
        end
    end
end

fprintf('Number of Screened Out Low-Viability Response Curves is : %d.\n', length(selected_LowViability_pubchemSID));
m = length(averaged_viability_data_with_SID_exclude_concentration(:,1)) - length(selected_LowViability_pubchemSID);
fprintf('Number of Remaining Response Curves is : %d.\n', m);
remaining_pubchemID_after_viability = setdiff(averaged_viability_data_with_SID_exclude_concentration(:,1),selected_LowViability_pubchemSID);
fprintf('Number of Remaining Active Chemicals after Viability is : %d.\n', length(remaining_pubchemID_after_viability));



% Identify chemicals active in luciferase inhibition assay (averaged intereference > -10)
selected_highIntereference_pubchemSID = [];
for i = 1:1:length(averaged_interference_data(:,1))
  obs = averaged_interference_data(i,2:length(averaged_interference_data(1,:)));

  if (min(obs) > interference_boundary)
  
    selected_highIntereference_pubchemSID = [selected_highIntereference_pubchemSID, averaged_interference_data(i,1)];
  end
end

fprintf('Number of Screened Out High-Interference Response Curves is : %d.\n', length(selected_highIntereference_pubchemSID));
m = length(averaged_interference_data(:,1)) - length(selected_highIntereference_pubchemSID);
fprintf('Number of Remaining Response Curves is : %d.\n', m);
remaining_pubchemID_after_interference = setdiff(averaged_interference_data(:,1),selected_highIntereference_pubchemSID);
fprintf('Number of Remaining Active Chemicals after Interference is : %d.\n', length(remaining_pubchemID_after_interference));



% Identify chemicals active in autofluoresecence assay (averaged autofluorescence < 10)
selected_lowAutofluorescence_pubchemSID = [];
for i = 1:1:length(averaged_autofluorescence_data(:,1))
  obs = averaged_autofluorescence_data(i,2:length(averaged_autofluorescence_data(1,:)));

  if (max(obs) < autofluorescence_boundary)
  
    selected_lowAutofluorescence_pubchemSID = [selected_lowAutofluorescence_pubchemSID, averaged_autofluorescence_data(i,1)];
  end
end

fprintf('Number of Screened Out Low-Autofluorescence Response Curves is : %d.\n', length(selected_lowAutofluorescence_pubchemSID));
m = length(averaged_autofluorescence_data(:,1)) - length(selected_lowAutofluorescence_pubchemSID);
fprintf('Number of Remaining Response Curves is : %d.\n', m);
remaining_pubchemID_after_autofluorescence = setdiff(averaged_autofluorescence_data(:,1),selected_lowAutofluorescence_pubchemSID);
fprintf('Number of Remaining Active Chemicals after autofluorescence is : %d.\n', length(remaining_pubchemID_after_autofluorescence));

% number of chemicals with both interference and autofluoresence
length(intersect(remaining_pubchemID_after_interference, remaining_pubchemID_after_autofluorescence))


%% Plot compounds active in luciferase inhibition and autofluoresecence assays

% luciferase inhibition
for i = 1:1:length(remaining_pubchemID_after_interference)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        remaining_pubchemID_after_interference(i));

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    plot(log10(conc),single_interference_data,'-o','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    hold on
end

% autofluoresence
for i = 1:1:length(remaining_pubchemID_after_autofluorescence)

    index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
        remaining_pubchemID_after_autofluorescence(i));

    conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
    plot(log10(conc),single_autofluorescence_data,'-o','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    hold on
end

for i = 1:1:length(selected_lowAutofluorescence_pubchemSID)

    index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
        selected_lowAutofluorescence_pubchemSID(i));

    conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
    plot(log10(conc),single_autofluorescence_data,'-o','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    hold on
end


%% (D) False-positive NMDR exclusion:(2) Filter out U-shaped CRCs for compounds that are active in cytotoxicity, luciferase inhibition, autofluoresecence assays 

flat_pubchem_ID_set = [];
increasing_pubchem_ID_set = [];
decreasing_pubchem_ID_set = [];
U_pubchem_ID_set = [];
Bell_pubchem_ID_set = [];
S_pubchem_ID_set = [];

for i = 1:1:length(Flat_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Flat_index(i),1);
    flat_pubchem_ID_set = [flat_pubchem_ID_set,pubchem_ID];  
end

for i = 1:1:length(Increasing_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Increasing_index(i),1);
    increasing_pubchem_ID_set = [increasing_pubchem_ID_set,pubchem_ID];  
end

for i = 1:1:length(Decreasing_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Decreasing_index(i),1);
    decreasing_pubchem_ID_set = [decreasing_pubchem_ID_set,pubchem_ID];  
end

for i = 1:1:length(U_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(U_index(i),1);
    U_pubchem_ID_set = [U_pubchem_ID_set,pubchem_ID];  
end

for i = 1:1:length(Bell_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(Bell_index(i),1);
    Bell_pubchem_ID_set = [Bell_pubchem_ID_set,pubchem_ID];  
end

for i = 1:1:length(S_index)
    pubchem_ID = averaged_response_data_exclude_low_variance_outstanding_peak(S_index(i),1);
    S_pubchem_ID_set = [S_pubchem_ID_set,pubchem_ID];  
end



% filter out U based on cytotoxicity
U_pubchemID_with_viability = [];
U_pubchemID_with_viability_correlation = [];
U_pubchemID_with_viability_correlation_drop = [];

for i = 1:1:length(U_pubchem_ID_set)
    for j = 1:1:length(remaining_pubchemID_after_viability)
     
        if(U_pubchem_ID_set(i) == remaining_pubchemID_after_viability(j))
            
            index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
                U_pubchem_ID_set(i));
            index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
                U_pubchem_ID_set(i));
            U_pubchemID_with_viability = [U_pubchemID_with_viability, U_pubchem_ID_set(i)];
            
            single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
            minimal_index = find (single_response_data == min(single_response_data));
            left_side_response_curve = single_response_data(1:minimal_index);
            
            single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
            left_side_viability_curve = single_viability_data(1:minimal_index);
            
            if (corr(left_side_response_curve',left_side_viability_curve') > viability_correlation_threshold)
                U_pubchemID_with_viability_correlation = [U_pubchemID_with_viability_correlation, U_pubchem_ID_set(i)];
                if((max(left_side_viability_curve)-min(left_side_viability_curve)) > viability_vertical_drop_threshold)
                    U_pubchemID_with_viability_correlation_drop = [U_pubchemID_with_viability_correlation_drop, U_pubchem_ID_set(i)];
                end
            end
            
        end
    end
end
U_pubchemID_without_viability = setdiff(U_pubchem_ID_set,U_pubchemID_with_viability)
length(U_pubchemID_without_viability) + length(U_pubchemID_with_viability)

for i = 1:1:length(U_pubchemID_with_viability)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability(i));

    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    left_side_viability_curve = single_viability_data(1:minimal_index);
    
    corr(left_side_response_curve',left_side_viability_curve')
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_pubchemID_with_viability(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_pubchemID_with_viability(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_pubchemID_with_viability(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(U_pubchemID_with_viability_correlation)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability_correlation(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability_correlation(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    left_side_viability_curve = single_viability_data(1:minimal_index);
    
    corr(left_side_response_curve',left_side_viability_curve')
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_pubchemID_with_viability_correlation(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_pubchemID_with_viability_correlation(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_pubchemID_with_viability_correlation(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(U_pubchemID_with_viability_correlation_drop)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability_correlation_drop(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_viability_correlation_drop(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    left_side_viability_curve = single_viability_data(1:minimal_index);
    
    corr(left_side_response_curve',left_side_viability_curve')
    (max(left_side_viability_curve)-min(left_side_viability_curve))

    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_pubchemID_with_viability_correlation_drop(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_pubchemID_with_viability_correlation_drop(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_pubchemID_with_viability_correlation_drop(i))],'Location','southwest');
    legend('boxoff')
end

% Plot remaining U curves after viability
final_U_curve_pubchemID_after_viability = setdiff(U_pubchem_ID_set,U_pubchemID_with_viability_correlation_drop);
length(final_U_curve_pubchemID_after_viability)
for i = 1:1:length(final_U_curve_pubchemID_after_viability)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID_after_viability(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID_after_viability(i));

    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    left_side_viability_curve = single_viability_data(1:minimal_index);
    
    corViability = corr(left_side_response_curve',left_side_viability_curve')
    dropViability = (max(left_side_viability_curve)-min(left_side_viability_curve))

    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_U_curve_pubchemID_after_viability(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == final_U_curve_pubchemID_after_viability(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)


    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_U_curve_pubchemID_after_viability(i))],'Location','southwest');
    text(-9.5, -10, num2str(corViability), 'Color', 'r', 'FontSize', 20);
    text(-9.5, -30, num2str(dropViability), 'Color', 'r', 'FontSize', 20);
    legend('boxoff')
end




% filter out U based on luciferase inhibition
U_pubchemID_with_interference = [];
U_pubchemID_with_interference_correlation = [];
U_pubchemID_with_interference_correlation_drop = [];
for i = 1:1:length(U_pubchem_ID_set)
    for j = 1:1:length(remaining_pubchemID_after_interference)
     
        if(U_pubchem_ID_set(i) == remaining_pubchemID_after_interference(j))
            U_pubchemID_with_interference = [U_pubchemID_with_interference, U_pubchem_ID_set(i)];  
             
            index_in_interference = find (averaged_interference_data(:,1) == ...
                U_pubchem_ID_set(i));
            index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
                U_pubchem_ID_set(i));
            
            single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
            minimal_index = find (single_response_data == min(single_response_data));
            left_side_response_curve = single_response_data(1:minimal_index);
            
            single_interference_data = averaged_interference_data(index_in_interference,2:end);
            left_side_interference_curve = single_interference_data(1:minimal_index);

            if (corr(left_side_response_curve',left_side_interference_curve') > interference_correlation_threshold)
                U_pubchemID_with_interference_correlation = [U_pubchemID_with_interference_correlation, U_pubchem_ID_set(i)];
                if((max(left_side_interference_curve)-min(left_side_interference_curve)) > interference_vertical_drop_threshold)
                    U_pubchemID_with_interference_correlation_drop = [U_pubchemID_with_interference_correlation_drop, U_pubchem_ID_set(i)];
                end
            end

        end  
    end
end

U_pubchemID_without_interference = setdiff(U_pubchem_ID_set,U_pubchemID_with_interference)
length(U_pubchemID_without_interference) + length(U_pubchemID_with_interference)
for i = 1:1:length(U_pubchemID_with_interference_correlation_drop)

    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        U_pubchemID_with_interference_correlation_drop(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_pubchemID_with_interference_correlation_drop(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    left_side_interference_curve = single_interference_data(1:minimal_index);
    
    corr(left_side_response_curve',left_side_interference_curve')
    max(left_side_interference_curve)-min(left_side_interference_curve)
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_pubchemID_with_interference_correlation_drop(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == U_pubchemID_with_interference_correlation_drop(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_pubchemID_with_interference_correlation_drop(i))],'Location','southwest');
    legend('boxoff')
end

% Remaining U curves after interference
final_U_curve_pubchemID_after_interference = setdiff(U_pubchem_ID_set,U_pubchemID_with_interference_correlation_drop);
length(final_U_curve_pubchemID_after_interference)

% Plot remaining U curves after interference without data
final_U_curve_pubchemID_after_interference_withoutdata = setdiff(final_U_curve_pubchemID_after_interference, averaged_interference_data(:,1));
length(final_U_curve_pubchemID_after_interference_withoutdata)
for i = 1:1:length(final_U_curve_pubchemID_after_interference_withoutdata)

    figure(i)

    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID_after_interference_withoutdata(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_U_curve_pubchemID_after_interference_withoutdata(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_U_curve_pubchemID_after_interference_withoutdata(i))],'Location','southwest');
    legend('boxoff')

end

% Plot remaining U curves after interference with data
final_U_curve_pubchemID_after_interference_withdata = setdiff(final_U_curve_pubchemID_after_interference,final_U_curve_pubchemID_after_interference_withoutdata);
length(final_U_curve_pubchemID_after_interference_withdata)
for i = 1:1:length(final_U_curve_pubchemID_after_interference_withdata)

    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        final_U_curve_pubchemID_after_interference_withdata(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID_after_interference_withdata(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    minimal_index = find (single_response_data == min(single_response_data));
    left_side_response_curve = single_response_data(1:minimal_index);
    
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    left_side_interference_curve = single_interference_data(1:minimal_index);
    
    corr(left_side_response_curve',left_side_interference_curve')
    max(left_side_interference_curve)-min(left_side_interference_curve)
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_U_curve_pubchemID_after_interference_withdata(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == final_U_curve_pubchemID_after_interference_withdata(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_U_curve_pubchemID_after_interference_withdata(i))],'Location','southwest');
    legend('boxoff')

end





% filter out U based on autofluoresecence
if (hasAutofluoresenceData == 1)
    U_pubchemID_with_autofluoresence = [];
    
    for i = 1:1:length(U_pubchem_ID_set)
        for j = 1:1:length(remaining_pubchemID_after_autofluorescence)
            
            if(U_pubchem_ID_set(i) == remaining_pubchemID_after_autofluorescence(j))
                U_pubchemID_with_autofluoresence = [U_pubchemID_with_autofluoresence, U_pubchem_ID_set(i)];
            end
        end
    end
    
    U_pubchemID_without_autofluoresence = setdiff(U_pubchem_ID_set,U_pubchemID_with_autofluoresence)
    length(U_pubchemID_without_autofluoresence) + length(U_pubchemID_with_autofluoresence)
    
    for i = 1:1:length(U_pubchemID_with_autofluoresence)
        
        figure(i)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            U_pubchemID_with_autofluoresence(i));
        index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
            U_pubchemID_with_autofluoresence(i));
        
        chemical_index_in_full_response_data = find(full_response_data(:,2) == U_pubchemID_with_autofluoresence(i));
        full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
        standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
        
        chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == U_pubchemID_with_autofluoresence(i));
        full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
        standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));
        
        conc = averaged_concentration_response(index_in_response,2:end);
        single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
        errorbar(log10(conc),single_response_data,...
            standard_deviation_for_response_data,'-ko','LineWidth',2)
        hold on
        
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-ro','LineWidth',2)
        %xlim([-4 3])
        %ylim([-100 100])
        %yticks([-100:20:100])
        xlabel('Dose')
        ylabel('Response')
        set(gca,'fontsize',18)
        box off
        grid on
        legend(['ID=' num2str(U_pubchemID_with_autofluoresence(i))],'Location','southwest');
        legend('boxoff')
    end
    
    % Visual inspection on U curve with autofluoresence
    selectedIndex = [4]; % index 1,2,3 was selected after visual inspection
    Removed_U_with_autofluoresence = [];
    for i = 1:1:length(selectedIndex)
        Removed_U_with_autofluoresence = [Removed_U_with_autofluoresence, U_pubchemID_with_autofluoresence(selectedIndex(i))];
    end
    
    % Plot remaining U curves after autofluoresence
    final_U_curve_pubchemID_after_autofluoresence = setdiff(U_pubchem_ID_set,Removed_U_with_autofluoresence);
    length(final_U_curve_pubchemID_after_autofluoresence)
    
    for i = 1:1:length(final_U_curve_pubchemID_after_autofluoresence)
        
        figure(i)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_U_curve_pubchemID_after_autofluoresence(i));
        index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
            final_U_curve_pubchemID_after_autofluoresence(i));
        
        single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        
        
        chemical_index_in_full_response_data = find(full_response_data(:,2) == final_U_curve_pubchemID_after_autofluoresence(i));
        full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
        standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
        
        
        chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == final_U_curve_pubchemID_after_autofluoresence(i));
        full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
        standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));
        
        
        conc = averaged_concentration_response(index_in_response,2:end);
        single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
        errorbar(log10(conc),single_response_data,...
            standard_deviation_for_response_data,'-ko','LineWidth',2)
        hold on
        
        
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-bo','LineWidth',2)
        
        
        
        %xlim([-4 3])
        %ylim([-100 100])
        %yticks([-100:20:100])
        xlabel('Dose')
        ylabel('Response')
        set(gca,'fontsize',18)
        box off
        grid on
        legend(['ID=' num2str(final_U_curve_pubchemID_after_autofluoresence(i))],'Location','southwest');
        legend('boxoff')
    end
end

    




% Plot final remainig U-shaped CRCs
if (hasAutofluoresenceData == 1)
    final_U_curve_pubchemID_temp = final_U_curve_pubchemID_after_autofluoresence;
end

if (hasIntereferenceData == 1)
    final_U_curve_pubchemID_temp = final_U_curve_pubchemID_after_interference;
end

final_U_curve_pubchemID = intersect(final_U_curve_pubchemID_temp,final_U_curve_pubchemID_after_viability);
length(final_U_curve_pubchemID)
index_in_interference = [];
index_in_autofluorescence = [];

for i = 1:1:length(final_U_curve_pubchemID)

    figure(i)
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_U_curve_pubchemID(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            final_U_curve_pubchemID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_U_curve_pubchemID(i));
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_U_curve_pubchemID(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == final_U_curve_pubchemID(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == final_U_curve_pubchemID(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == final_U_curve_pubchemID(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data/sqrt(length(chemical_index_in_full_response_data)),'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data/sqrt(length(chemical_index_in_full_viability_data)),'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data/sqrt(length(chemical_index_in_full_interference_data)),'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data/sqrt(length(chemical_index_in_full_autofluoresence_data)),'-bo','LineWidth',2)
    end
    
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_U_curve_pubchemID(i))],'Location','southwest');
    legend('boxoff')
end

% Plot U curves filtered out due to interference and/or autofluoresence
U_curve_filtered_out_dueTo_interference_autofluoresence = setdiff(final_U_curve_pubchemID_after_viability,final_U_curve_pubchemID);
for i = 1:1:length(U_curve_filtered_out_dueTo_interference_autofluoresence)

    figure(i)
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_interference_autofluoresence(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_interference_autofluoresence(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_curve_filtered_out_dueTo_interference_autofluoresence(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_curve_filtered_out_dueTo_interference_autofluoresence(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == U_curve_filtered_out_dueTo_interference_autofluoresence(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == U_curve_filtered_out_dueTo_interference_autofluoresence(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data,'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-bo','LineWidth',2)
    end
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_curve_filtered_out_dueTo_interference_autofluoresence(i))],'Location','southwest');
    legend('boxoff')
end

% Plot U curves filtered out due to viability and/or autofluoresence
U_curve_filtered_out_dueTo_viability_autofluoresence = setdiff(final_Bell_curve_pubchemID_after_interference,final_U_curve_pubchemID);
for i = 1:1:length(U_curve_filtered_out_dueTo_viability_autofluoresence)

    figure(i)
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_viability_autofluoresence(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_viability_autofluoresence(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_curve_filtered_out_dueTo_viability_autofluoresence(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_curve_filtered_out_dueTo_viability_autofluoresence(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == U_curve_filtered_out_dueTo_viability_autofluoresence(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == U_curve_filtered_out_dueTo_viability_autofluoresence(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data,'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-bo','LineWidth',2)
    end
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_curve_filtered_out_dueTo_viability_autofluoresence(i))],'Location','southwest');
    legend('boxoff')
end

% Plot U curves filtered out due to viability and/or interference
U_curve_filtered_out_dueTo_viability_intereference = setdiff(final_U_curve_pubchemID_after_autofluoresence,final_U_curve_pubchemID);
for i = 1:1:length(U_curve_filtered_out_dueTo_viability_intereference)

    figure(i)
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_viability_intereference(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        U_curve_filtered_out_dueTo_viability_intereference(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == U_curve_filtered_out_dueTo_viability_intereference(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == U_curve_filtered_out_dueTo_viability_intereference(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == U_curve_filtered_out_dueTo_viability_intereference(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == U_curve_filtered_out_dueTo_viability_intereference(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data,'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-bo','LineWidth',2)
    end
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(U_curve_filtered_out_dueTo_viability_intereference(i))],'Location','southwest');
    legend('boxoff')
end


%% (D) False-positive NMDR exclusion:(3) Filter out Bell-shaped CRCs for compounds that are active in cytotoxicity, luciferase inhibition, autofluoresecence assays

% Filter out Bell based on cytotoxicity
Bell_pubchemID_with_viability = [];
Bell_pubchemID_with_viability_correlation = [];
Bell_pubchemID_with_viability_correlation_drop = [];

for i = 1:1:length(Bell_pubchem_ID_set)
    for j = 1:1:length(remaining_pubchemID_after_viability)

        if(Bell_pubchem_ID_set(i) == remaining_pubchemID_after_viability(j))

            index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
                Bell_pubchem_ID_set(i));
            index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
                Bell_pubchem_ID_set(i));
            Bell_pubchemID_with_viability = [Bell_pubchemID_with_viability, Bell_pubchem_ID_set(i)];
            
            single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
            maximal_index = find (single_response_data == max(single_response_data));
            right_side_response_curve = single_response_data(maximal_index:end);
            
            single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
            right_side_viability_curve = single_viability_data(maximal_index:end);
            
            if (corr(right_side_response_curve',right_side_viability_curve') > viability_correlation_threshold)
                Bell_pubchemID_with_viability_correlation = [Bell_pubchemID_with_viability_correlation, Bell_pubchem_ID_set(i)];
                if((max(right_side_viability_curve)-min(right_side_viability_curve)) > viability_vertical_drop_threshold)
                    Bell_pubchemID_with_viability_correlation_drop = [Bell_pubchemID_with_viability_correlation_drop, Bell_pubchem_ID_set(i)];
                end
            end
            
        end
    end
end
Bell_pubchemID_without_viability = setdiff(Bell_pubchem_ID_set,Bell_pubchemID_with_viability)
length(Bell_pubchemID_without_viability) + length(Bell_pubchemID_with_viability)

for i = 1:1:length(Bell_pubchemID_with_viability)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    right_side_viability_curve = single_viability_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_viability_curve')
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_viability(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == Bell_pubchemID_with_viability(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_viability(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(Bell_pubchemID_with_viability_correlation)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability_correlation(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability_correlation(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    right_side_viability_curve = single_viability_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_viability_curve')
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_viability_correlation(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == Bell_pubchemID_with_viability_correlation(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_viability_correlation(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(Bell_pubchemID_with_viability_correlation_drop)

    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability_correlation_drop(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_viability_correlation_drop(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    right_side_viability_curve = single_viability_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_viability_curve')
    (max(right_side_viability_curve)-min(right_side_viability_curve))

    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_viability_correlation_drop(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == Bell_pubchemID_with_viability_correlation_drop(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_viability_correlation_drop(i))],'Location','southwest');
    legend('boxoff')
end

% Plot remaining Bell curves after viability

final_Bell_curve_pubchemID_after_viability = setdiff(Bell_pubchem_ID_set,Bell_pubchemID_with_viability_correlation_drop);
length(final_Bell_curve_pubchemID_after_viability)
for i = 1:1:length(final_Bell_curve_pubchemID_after_viability)
    
    figure(i)
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID_after_viability(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID_after_viability(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    right_side_viability_curve = single_viability_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_viability_curve')
    max(right_side_viability_curve)-min(right_side_viability_curve)

    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_Bell_curve_pubchemID_after_viability(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == final_Bell_curve_pubchemID_after_viability(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_Bell_curve_pubchemID_after_viability(i))],'Location','southwest');
    legend('boxoff')
end





% Filter out Bell based on luciferase inhibition
Bell_pubchemID_with_interference = [];
Bell_pubchemID_with_interference_correlation = [];
Bell_pubchemID_with_interference_correlation_drop = [];

for i = 1:1:length(Bell_pubchem_ID_set)
    for j = 1:1:length(remaining_pubchemID_after_interference)
     
        if(Bell_pubchem_ID_set(i) == remaining_pubchemID_after_interference(j))
            Bell_pubchemID_with_interference = [Bell_pubchemID_with_interference, Bell_pubchem_ID_set(i)];  
             
            index_in_interference = find (averaged_interference_data(:,1) == ...
                Bell_pubchem_ID_set(i));
            index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
                Bell_pubchem_ID_set(i));
            
            single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
            maximal_index = find (single_response_data == max(single_response_data));
            right_side_response_curve = single_response_data(maximal_index:end);
            
            single_interference_data = averaged_interference_data(index_in_interference,2:end);
            right_side_interference_curve = single_interference_data(maximal_index:end);

            if (corr(right_side_response_curve',right_side_interference_curve') > interference_correlation_threshold)
                Bell_pubchemID_with_interference_correlation = [Bell_pubchemID_with_interference_correlation, Bell_pubchem_ID_set(i)];
                if((max(right_side_interference_curve)-min(right_side_interference_curve)) > interference_vertical_drop_threshold)
                    Bell_pubchemID_with_interference_correlation_drop = [Bell_pubchemID_with_interference_correlation_drop, Bell_pubchem_ID_set(i)];
                end
            end

        end  
    end
end

Bell_pubchemID_without_interference = setdiff(Bell_pubchem_ID_set,Bell_pubchemID_with_interference)
length(Bell_pubchemID_without_interference) + length(Bell_pubchemID_with_interference)

for i = 1:1:length(Bell_pubchemID_with_interference)
    
    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        Bell_pubchemID_with_interference(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_interference(i));
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_interference(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == Bell_pubchemID_with_interference(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_interference(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(Bell_pubchemID_with_interference_correlation)

    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        Bell_pubchemID_with_interference_correlation(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_interference_correlation(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    right_side_interference_curve = single_interference_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_interference_curve')
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_interference_correlation(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == Bell_pubchemID_with_interference_correlation(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_interference_correlation(i))],'Location','southwest');
    legend('boxoff')
end

for i = 1:1:length(Bell_pubchemID_with_interference_correlation_drop)

    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        Bell_pubchemID_with_interference_correlation_drop(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        Bell_pubchemID_with_interference_correlation_drop(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    right_side_interference_curve = single_interference_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_interference_curve')
    max(right_side_interference_curve)-min(right_side_interference_curve)
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_interference_correlation_drop(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == Bell_pubchemID_with_interference_correlation_drop(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(Bell_pubchemID_with_interference_correlation_drop(i))],'Location','southwest');
    legend('boxoff')
end

% Remaining Bell curves after interference
final_Bell_curve_pubchemID_after_interference = setdiff(Bell_pubchem_ID_set,Bell_pubchemID_with_interference_correlation_drop);
length(final_Bell_curve_pubchemID_after_interference)

% Plot remaining Bell curves after interference without data
final_Bell_curve_pubchemID_after_interference_withoutdata = setdiff(final_Bell_curve_pubchemID_after_interference, averaged_interference_data(:,1));
length(final_Bell_curve_pubchemID_after_interference_withoutdata)
for i = 1:1:length(final_Bell_curve_pubchemID_after_interference_withoutdata)

    figure(i)

    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID_after_interference(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_Bell_curve_pubchemID_after_interference(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_Bell_curve_pubchemID_after_interference(i))],'Location','southwest');
    legend('boxoff')

end

% Plot remaining Bell curves after interference with data
final_Bell_curve_pubchemID_after_interference_withdata = setdiff(final_Bell_curve_pubchemID_after_interference,final_Bell_curve_pubchemID_after_interference_withoutdata);
length(final_Bell_curve_pubchemID_after_interference_withdata)
for i = 1:1:length(final_Bell_curve_pubchemID_after_interference_withdata)

    figure(i)
    index_in_interference = find (averaged_interference_data(:,1) == ...
        final_Bell_curve_pubchemID_after_interference_withdata(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID_after_interference_withdata(i));
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    maximal_index = find (single_response_data == max(single_response_data));
    right_side_response_curve = single_response_data(maximal_index:end);
    
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    right_side_interference_curve = single_interference_data(maximal_index:end);
    
    corr(right_side_response_curve',right_side_interference_curve')
    max(right_side_interference_curve)-min(right_side_interference_curve)
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_Bell_curve_pubchemID_after_interference_withdata(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == final_Bell_curve_pubchemID_after_interference_withdata(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data,'-ko','LineWidth',2)
    hold on

    conc = averaged_concentration_for_interference(index_in_interference,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    errorbar(log10(conc),single_interference_data,...
        standard_deviation_for_interference_data,'-ro','LineWidth',2)
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_Bell_curve_pubchemID_after_interference_withdata(i))],'Location','southwest');
    legend('boxoff')

end





% Filter out Bell based on autofluoresecence 
if (hasAutofluoresenceData == 1)
    Bell_pubchemID_with_autofluoresence = [];
    
    for i = 1:1:length(Bell_pubchem_ID_set)
        for j = 1:1:length(remaining_pubchemID_after_autofluorescence)
            
            if(Bell_pubchem_ID_set(i) == remaining_pubchemID_after_autofluorescence(j))
                Bell_pubchemID_with_autofluoresence = [Bell_pubchemID_with_autofluoresence, Bell_pubchem_ID_set(i)];
            end
        end
    end
    
    Bell_pubchemID_without_autofluoresence = setdiff(Bell_pubchem_ID_set,Bell_pubchemID_with_autofluoresence)
    length(Bell_pubchemID_without_autofluoresence) + length(Bell_pubchemID_with_autofluoresence)
    
    for i = 1:1:length(Bell_pubchemID_with_autofluoresence)
        
        figure(i)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            Bell_pubchemID_with_autofluoresence(i));
        index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
            Bell_pubchemID_with_autofluoresence(i));
        
        chemical_index_in_full_response_data = find(full_response_data(:,2) == Bell_pubchemID_with_autofluoresence(i));
        full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
        standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
        
        chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == Bell_pubchemID_with_autofluoresence(i));
        full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
        standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));
        
        conc = averaged_concentration_response(index_in_response,2:end);
        single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
        errorbar(log10(conc),single_response_data,...
            standard_deviation_for_response_data,'-ko','LineWidth',2)
        hold on
        
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data,'-ro','LineWidth',2)
        %xlim([-4 3])
        %ylim([-100 100])
        %yticks([-100:20:100])
        xlabel('Dose')
        ylabel('Response')
        set(gca,'fontsize',18)
        box off
        grid on
        legend(['ID=' num2str(Bell_pubchemID_with_autofluoresence(i))],'Location','southwest');
        legend('boxoff')
    end
    
    % Visual inspection on U curve with autofluoresence
    selectedIndex = [1,3,4,6]; % index 1,2,3 was selected after visual inspection
    Removed_Bell_with_autofluoresence = [];
    for i = 1:1:length(selectedIndex)
        Removed_Bell_with_autofluoresence = [Removed_Bell_with_autofluoresence, Bell_pubchemID_with_autofluoresence(selectedIndex(i))];
    end
    
    % Remaining Bell curves after autofluoresence
    final_Bell_curve_pubchemID_after_autofluoresence = setdiff(Bell_pubchem_ID_set,Removed_Bell_with_autofluoresence);
    length(final_Bell_curve_pubchemID_after_autofluoresence)
end






% Plot final remainig Bell-shaped CRCs
if (hasAutofluoresenceData == 1)
    final_Bell_curve_pubchemID_temp = final_Bell_curve_pubchemID_after_autofluoresence;
end

if (hasIntereferenceData == 1)
    final_Bell_curve_pubchemID_temp = final_Bell_curve_pubchemID_after_interference;
end

final_Bell_curve_pubchemID = intersect(final_Bell_curve_pubchemID_temp,final_Bell_curve_pubchemID_after_viability);
length(final_Bell_curve_pubchemID)
index_in_interference = [];
index_in_autofluorescence = [];


for i = 1:1:length(final_Bell_curve_pubchemID)

    figure(i)

    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID(i));
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        final_Bell_curve_pubchemID(i));
   
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            final_Bell_curve_pubchemID(i));
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);

    chemical_index_in_full_response_data = find(full_response_data(:,2) == final_Bell_curve_pubchemID(i));
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == final_Bell_curve_pubchemID(i));
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));

    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == final_Bell_curve_pubchemID(i));
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == final_Bell_curve_pubchemID(i));
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));

    conc = averaged_concentration_response(index_in_response,2:end);
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    errorbar(log10(conc),single_response_data,...
        standard_deviation_for_response_data/sqrt(length(chemical_index_in_full_response_data)),'-ko','LineWidth',2)
    hold on
    
    conc = averaged_concentration_viability(index_in_viability,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    errorbar(log10(conc),single_viability_data,...
        standard_deviation_for_viability_data/sqrt(length(chemical_index_in_full_viability_data)),'-ro','LineWidth',2)
    
    hold on
    
    if (~isempty(index_in_interference))
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data/sqrt(length(chemical_index_in_full_interference_data)),'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence))
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data/sqrt(length(chemical_index_in_full_autofluoresence_data)),'-bo','LineWidth',2)
    end
    
    
    %xlim([-4 3])
    %ylim([-100 100])
    %yticks([-100:20:100])
    xlabel('Dose')
    ylabel('Response')
    set(gca,'fontsize',18)
    box off
    grid on
    legend(['ID=' num2str(final_Bell_curve_pubchemID(i))],'Location','southwest');
    legend('boxoff')
end


%% Scatter plots of log10 peak-activity concentration vs. magnitude of identified NMDR CRCs

% X axis represents the concentration where the peak or nadir point appears (log nM)
% Y axis represents the magnitude of Bell or U-shaped CRCs 

volcano_matrix_U = [];
volcano_matrix_Bell = [];
volcano_matrix_increasing = [];

for i = 1:1:length(final_U_curve_pubchemID)
    % process response data

    index_in_response = find (averaged_response_data_exclude_low_variance_outstanding_peak(:,1) == final_U_curve_pubchemID(i));
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(index_in_response,:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(index_in_response,:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(index_in_response,:);

     
%     single_response_exclude_NA = single_response_data;
%     single_response_exclude_NA = single_response_exclude_NA(~isnan(single_response_exclude_NA)); 
%     single_response_MA_1_exclude_NA = single_response_MA_1;
%     single_response_MA_1_exclude_NA = single_response_MA_1_exclude_NA(~isnan(single_response_MA_1_exclude_NA)); 
%     single_response_MA_2_exclude_NA = single_response_MA_2;
%     single_response_MA_2_exclude_NA = single_response_MA_2_exclude_NA(~isnan(single_response_MA_2_exclude_NA)); 
 
    % When calculate the curvature, use the data after moving average;
    % because the U index is obtained based on those data;
    % It is possible that the imputated data is identified as a U shape but
    % the original data is not
    min_index = find(single_response_data == min(single_response_data));
    first_index = 1;
    last_index = length(single_response_data);
    curvature_degree = min([abs(single_response_data(min_index)-single_response_data(first_index)),...
        abs(single_response_data(min_index)-single_response_data(last_index))]);
    % the concentration where inflexion point appears is based on the
    % original data; not the data excluding NAN values

    conc = averaged_concentration_response(index_in_response,2:end);
    conc = log10(conc);
    conc_index = find(single_response_data == min(single_response_data));
    conc_at_inflexion = conc(conc_index);
    volcano_matrix_U(i,1) = conc_at_inflexion;
    volcano_matrix_U(i,2) = curvature_degree;
    volcano_matrix_U(i,3) = final_U_curve_pubchemID(i);
end

for i = 1:1:length(final_Bell_curve_pubchemID)
    % process response data

    index_in_response = find (averaged_response_data_exclude_low_variance_outstanding_peak(:,1) == final_Bell_curve_pubchemID(i));
    single_response_data = AveragedResponseDataExcludeLowVarianceOutstandingPeakPubchemSID(index_in_response,:);
    single_response_MA_1 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_1(index_in_response,:);
    single_response_MA_2 = AveragedResponseDataExcludeLowVariancePeakPubchemSID_MA_2(index_in_response,:);

     
%     single_response_exclude_NA = single_response_data;
%     single_response_exclude_NA = single_response_exclude_NA(~isnan(single_response_exclude_NA)); 
%     single_response_MA_1_exclude_NA = single_response_MA_1;
%     single_response_MA_1_exclude_NA = single_response_MA_1_exclude_NA(~isnan(single_response_MA_1_exclude_NA)); 
%     single_response_MA_2_exclude_NA = single_response_MA_2;
%     single_response_MA_2_exclude_NA = single_response_MA_2_exclude_NA(~isnan(single_response_MA_2_exclude_NA)); 
 
    % When calculate the curvature, use the data after moving average;
    % because the U index is obtained based on those data;
    % It is possible that the imputated data is identified as a U shape but
    % the original data is not
    max_index = find(single_response_data == max(single_response_data));
    first_index = 1;
    last_index = length(single_response_data);
    curvature_degree = min([abs(single_response_data(max_index)-single_response_data(first_index)),...
        abs(single_response_data(max_index)-single_response_data(last_index))]);
    % the concentration where inflexion point appears is based on the
    % original data; not the data excluding NAN values

    conc = averaged_concentration_response(index_in_response,2:end);
    conc = log10(conc);
    conc_index = find(single_response_data == max(single_response_data));
    conc_at_inflexion = conc(conc_index);
    volcano_matrix_Bell(i,1) = conc_at_inflexion;
    volcano_matrix_Bell(i,2) = curvature_degree;
    volcano_matrix_Bell(i,3) = final_Bell_curve_pubchemID(i);
end


% Scatter plot for U-shaped CRCs 
figure(100)
plot(volcano_matrix_U(:,1),volcano_matrix_U(:,2),'bo','LineWidth',2)
xlim([-10 -4])
xticks([-10:2:-4])
% xlabel('Curvature')
% ylabel('-log10 Concentration')
set(gca,'fontsize',18)
box off
grid on


% Scatter plot for Bell-shaped CRCs 
figure(110)
plot(volcano_matrix_Bell(:,1),volcano_matrix_Bell(:,2),'bo','LineWidth',2)
xlim([-10 -4])
xticks([-10:2:-4])
% xlabel('Curvature')
% ylabel('-log10 Concentration')
set(gca,'fontsize',18)
box off
grid on

