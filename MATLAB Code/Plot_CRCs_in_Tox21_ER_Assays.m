% MatLab code for "Identification of nonmonotonic concentration-responses
% in Tox21 high- throughput screening estrogen receptor assays" published
% in Toxicology and Applied Pharmcology 2022

% This MatLab code is used to display a compound's concentration-response curves (CRCs) in all eight Tox21 ER assays given its PubChem SID (substance identifier):
% Note: Assays 2 and 5 do not contain cytotoxicity screening data; the cytotoxicity data of assays 4 and 6 were used instead for assays 2 and 5 respectively because they used the same cell reporter systems and identical compound library. 
clear all
clc
close all

% Type in compound PubChem SID 
SID = 144214049; %Bisphenol A

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


%% Plot concentration-response curves, cytotoxicity curves, luciferase inhibition curves, autofluoresecence curves in all ER assays where available

SID_chemical_name = readcell('Tox21_13127_Compounds_Library.csv');

% find compound index based on given SID
for i = 2:1:length(SID_chemical_name(:,1))
    if(cell2mat(SID_chemical_name(i,2)) == SID)
        compound_index = i;
    end
end
% find chemical name based on the compound index
chemical_name = SID_chemical_name(compound_index,3);
chemical_name = string(chemical_name); % convert chemical_name to a string variable

index_in_autofluorescence = [];
index_in_interference  = [];
figure(100)

for j = 1:1:8
    
    if (j == 1)
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        averaged_concentration_response = averaged_concentration_response_er_luc_agonist_p4;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_agonist_p4;
        full_response_data = full_response_data_er_luc_agonist_p4;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_agonist_p4;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_agonist_p4;
        full_viability_data = full_viability_data_er_luc_agonist_p4;
        
    end
    
    
    if (j == 2)
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        averaged_concentration_response = averaged_concentration_response_er_luc_agonist_p2;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_agonist_p2;
        full_response_data = full_response_data_er_luc_agonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p1;
        full_viability_data = full_viability_data_er_luc_antagonist_p1;
        
    end
    
    if (j == 3)
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        averaged_concentration_response = averaged_concentration_response_er_luc_antagonist_p2;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_antagonist_p2;
        full_response_data = full_response_data_er_luc_antagonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p2;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p2;
        full_viability_data = full_viability_data_er_luc_antagonist_p2;
        
    end
    
    
    if (j == 4)
        hasIntereferenceData = 1;
        hasAutofluoresenceData = 0;
        averaged_concentration_response = averaged_concentration_response_er_luc_antagonist_p1;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_luc_antagonist_p1;
        full_response_data = full_response_data_er_luc_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_er_luc_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_luc_antagonist_p1;
        full_viability_data = full_viability_data_er_luc_antagonist_p1;
        
    end
    
    
    if (j == 5)
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        averaged_concentration_response = averaged_concentration_response_er_bla_agonist_p2;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_bla_agonist_p2;
        full_response_data = full_response_data_er_bla_agonist_p2;
        averaged_concentration_viability = averaged_concentration_viability_er_bla_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_bla_antagonist_p1;
        full_viability_data = full_viability_data_er_bla_antagonist_p1;
        
    end
    
    
    if (j == 6)
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        averaged_concentration_response = averaged_concentration_response_er_bla_antagonist_p1;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_er_bla_antagonist_p1;
        full_response_data = full_response_data_er_bla_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_er_bla_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_er_bla_antagonist_p1;
        full_viability_data = full_viability_data_er_bla_antagonist_p1;
        
    end
    
    
    if (j == 7)
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        averaged_concentration_response = averaged_concentration_response_erb_bla_p1;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_erb_bla_p1;
        full_response_data = full_response_data_erb_bla_p1;
        averaged_concentration_viability = averaged_concentration_viability_erb_bla_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_erb_bla_p1;
        full_viability_data = full_viability_data_erb_bla_p1;
        
    end
    
    
    if (j == 8)
        hasIntereferenceData = 0;
        hasAutofluoresenceData = 1;
        averaged_concentration_response = averaged_concentration_response_erb_bla_antagonist_p1;
        averaged_response_data_with_SID_exclude_concentration = averaged_response_data_erb_bla_antagonist_p1;
        full_response_data = full_response_data_erb_bla_antagonist_p1;
        averaged_concentration_viability = averaged_concentration_viability_erb_bla_antagonist_p1;
        averaged_viability_data_with_SID_exclude_concentration = averaged_viability_data_erb_bla_antagonist_p1;
        full_viability_data = full_viability_data_erb_bla_antagonist_p1;
        
    end
    
    
    index_in_viability = find (averaged_viability_data_with_SID_exclude_concentration(:,1) == ...
        SID);
    index_in_response = find (averaged_response_data_with_SID_exclude_concentration(:,1) == ...
        SID);
    
    if (hasIntereferenceData == 1)
        index_in_interference = find (averaged_interference_data(:,1) == ...
            SID);
    end
    
    if (hasAutofluoresenceData == 1)
        index_in_autofluorescence = find (averaged_autofluorescence_data(:,1) == ...
            SID);
    end
    
    single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
    single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
    single_interference_data = averaged_interference_data(index_in_interference,2:end);
    single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
    
    chemical_index_in_full_response_data = find(full_response_data(:,2) == SID);
    full_response_data_one_chemical = full_response_data(chemical_index_in_full_response_data,3:17);
    standard_deviation_for_response_data = sqrt(var(full_response_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_viability_data = find(full_viability_data(:,2) == SID);
    full_viability_data_one_chemical = full_viability_data(chemical_index_in_full_viability_data,3:17);
    standard_deviation_for_viability_data = sqrt(var(full_viability_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_interference_data = find(full_interference_data(:,2) == SID);
    full_interference_data_one_chemical = full_interference_data(chemical_index_in_full_interference_data,3:17);
    standard_deviation_for_interference_data = sqrt(var(full_interference_data_one_chemical, 'omitnan'));
    
    chemical_index_in_full_autofluoresence_data = find(full_autofluorescence_data(:,2) == SID);
    full_autofluorescence_data_one_chemical = full_autofluorescence_data(chemical_index_in_full_autofluoresence_data,3:17);
    standard_deviation_for_autofluorescence_data = sqrt(var(full_autofluorescence_data_one_chemical, 'omitnan'));
    
%title(chemical_name);
    subplot(4,2,j);

    if(~isempty(index_in_response))
        conc = averaged_concentration_response(index_in_response,2:end);
        single_response_data = averaged_response_data_with_SID_exclude_concentration(index_in_response,2:end);
        errorbar(log10(conc),single_response_data,...
            standard_deviation_for_response_data/sqrt(length(chemical_index_in_full_response_data)),'-ko','LineWidth',2)
        hold on
    end
    
    
    if(~isempty(index_in_viability))
        conc = averaged_concentration_viability(index_in_viability,2:end);
        single_viability_data = averaged_viability_data_with_SID_exclude_concentration(index_in_viability,2:end);
        errorbar(log10(conc),single_viability_data,...
            standard_deviation_for_viability_data/sqrt(length(chemical_index_in_full_viability_data)),'-ro','LineWidth',2)
        
        hold on
    end
    
    if (~isempty(index_in_interference)&& hasIntereferenceData == 1)
        conc = averaged_concentration_for_interference(index_in_interference,2:end);
        single_interference_data = averaged_interference_data(index_in_interference,2:end);
        errorbar(log10(conc),single_interference_data,...
            standard_deviation_for_interference_data/sqrt(length(chemical_index_in_full_interference_data)),'-go','LineWidth',2)
    end
    
    hold on
    
    if (~isempty(index_in_autofluorescence) && hasAutofluoresenceData == 1)
        conc = averaged_concentration_for_autofluorescence(index_in_autofluorescence,2:end);
        single_autofluorescence_data = averaged_autofluorescence_data(index_in_autofluorescence,2:end);
        errorbar(log10(conc),single_autofluorescence_data,...
            standard_deviation_for_autofluorescence_data/sqrt(length(chemical_index_in_full_autofluoresence_data)),'-bo','LineWidth',2)
    end
    
    if j==7 || j==8
        xlabel('Log10 Concentration (M)')
    end
    ylabel('Activity Level')
    set(gca,'fontsize',12)
    box off
    grid on
    pbaspect([3,1,1])
    
    sample_number = length(chemical_index_in_full_response_data);
    
    legend off
    if(hasIntereferenceData == 1 & j==1)
        legend ('Activity','Cytotoxicity','Luciferase inhibition', 'Location','northwest')
        legend('boxoff')
    end
    if(hasAutofluoresenceData == 1 & j==5)
        legend ('Activity','Cytotoxicity','Autofluoresecence', 'Location','northwest')
        legend('boxoff')
    end

    title(['Assay ' num2str(j)]);

      
end

sgtitle(strcat(chemical_name, ' (SID: ', num2str(SID), ')'))



