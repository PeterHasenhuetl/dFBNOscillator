% Analysis of dFBN connectome (hemibrain:v1.2.1 and FlyWire FAFB v783 datasets). 
% Code written by Peter Hasenhuetl.


clear all
save_cond = 1;
if save_cond == 1

    source_data_details.data_path = []; %Add path as character array
    source_data_details.file_name_hemibrain = 'connectomics_analysis_hemibrain.xlsx';
    source_data_details.file_name_flywire = 'connectomics_analysis_flywire.xlsx';
    source_data_details.file_name_MB = 'connectomics_analysis_MB.xlsx';

else

    source_data_details = [];

end

conn_table_dFB_hemibrain = readtable('All_dFB_to_dFB_connections.csv');

if save_cond == 1
    writetable(conn_table_dFB_hemibrain, ...
        [source_data_details.data_path, source_data_details.file_name_hemibrain],...
        'Sheet', 'individual_connections')
end

[conn_table_dFB_flywire, details_array_dFB_flywire] = get_FlyWire_conn_table(source_data_details);
[conn_table_MB_flywire, adj_matrix_MB, symmetry_and_density_measures_MB] = get_MB_conn_table(source_data_details);

color = get_color;
y_pos = 4.5;
close all
figure('Name','connectomics analysis','Color','white',...
    'Units','centimeters','Position',[10 12 18 18],'Resize','off')

[adj_matrix_hemibrain, symmetry_and_density_measures_hemibrain] = ...
    get_connectome_analysis_and_fig(conn_table_dFB_hemibrain, y_pos+7, color, source_data_details, "hemibrain");
[adj_matrix_flywire, symmetry_and_density_measures_flywire] = ...
    get_connectome_analysis_and_fig(conn_table_dFB_flywire, y_pos, color, source_data_details, "flywire");
get_fig_panel_reciprocity(2, y_pos-3.5, symmetry_and_density_measures_hemibrain,...
    symmetry_and_density_measures_flywire, symmetry_and_density_measures_MB, color)
get_fig_panel_num_symmetry(10, y_pos-3.5, symmetry_and_density_measures_hemibrain,...
    symmetry_and_density_measures_flywire, symmetry_and_density_measures_MB, color)

if isempty(source_data_details) == 0
    cd(source_data_details.data_path)
    set(gcf,'renderer','Painters')
    saveas(gcf,'connectomics_analysis_fig.pdf')
    cd([]) %Add path as character array
end

%% Custom function to analyse and plot individual connectome

function [adj_matrix_1, symmetry_and_density_measures] = get_connectome_analysis_and_fig(...
    conn_table, y_pos, color, source_data_details, connectome_cond)

[directed_table, undirected_table, ~, adj_matrix_1, adj_matrix_sorted, ~, ~,...
    idx_left, idx_right, a_left, a_right] = get_undirected_graph(conn_table, 'all',...
    source_data_details, connectome_cond);
symmetry_and_density_measures = get_symmetry_and_density_measures(adj_matrix_1, source_data_details, connectome_cond);
symmetry_and_density_measures.adj_matrix_sorted = adj_matrix_sorted;

% Examines the weights between cell types and within cell types
dFB_layer_pre = [];
dFB_layer_post = [];
[ipsi_weights_within, contra_weights_within, size_ipsi_within, size_contra_within] = ...
    get_conn_weights(conn_table, dFB_layer_pre, dFB_layer_post, "within", "no");
[ipsi_weights_between, contra_weights_between, size_ipsi_between, size_contra_between] = ...
    get_conn_weights(conn_table, dFB_layer_pre, dFB_layer_post, "between", "no");

% Examines the weights between layers and within layers
[ipsi_weights_6_to_7, contra_weights_6_to_7, ~, ~] = ...
    get_conn_weights(conn_table, '6', '7', "all", "yes");
[ipsi_weights_7_to_6, contra_weights_7_to_6, ~, ~] = ...
    get_conn_weights(conn_table, '7', '6', "all", "yes");
inter_layer_weights = [ipsi_weights_6_to_7; contra_weights_6_to_7; ipsi_weights_7_to_6; contra_weights_7_to_6];

[ipsi_weights_6_to_6, contra_weights_6_to_6, ~, ~] = ...
    get_conn_weights(conn_table, '6', '6', "all", "yes");
[ipsi_weights_7_to_7, contra_weights_7_to_7, ~, ~] = ...
    get_conn_weights(conn_table, '7', '7', "all", "yes");
intra_layer_weights = [ipsi_weights_6_to_6; contra_weights_6_to_6; ipsi_weights_7_to_7; contra_weights_7_to_7];

% Calls the plotting functions
get_fig_panel_adjacency_matrix_heatmap(2, y_pos, adj_matrix_sorted, ...
    idx_left, idx_right, a_left, a_right, color, connectome_cond)
get_fig_panel_log_log_dist_v2(directed_table, undirected_table, color, 10, y_pos+3.2, source_data_details, connectome_cond)
ylim_1 = [0,0.8];
major_y_ticks = 0.4;
get_fig_panel_weight_histogram(inter_layer_weights, intra_layer_weights, ...
    14, y_pos+3.2, color, ylim_1, major_y_ticks, size_ipsi_within(1), ...
    size_contra_within(1),'Fraction of connections', 'layer_comparison', source_data_details,"layers", connectome_cond)

ylim_1 = [0,0.75];    
major_y_ticks = 0.25;
get_fig_panel_weight_histogram(ipsi_weights_between, contra_weights_between, ...
    10, y_pos, color, ylim_1, major_y_ticks, size_ipsi_between(1), ...
    size_contra_between(1), {'Fraction of connections';'between different cell types'},...
    'different_celltype', source_data_details, "hemispheres", connectome_cond)

get_fig_panel_weight_histogram(ipsi_weights_within, contra_weights_within, ...
    14, y_pos, color, ylim_1, major_y_ticks, size_ipsi_within(1), ...
    size_contra_within(1), {'Fraction of connections';'within same cell type'},...
    'same_celltype', source_data_details, "hemispheres", connectome_cond)

end

%% Custom analysis functions

function [directed_graph_edges, undirected_graph_edges, conn_table, ...
    adj_matrix, adj_matrix_sorted, curr_name, curr_ID,...
    idx_left, idx_right, a_left, a_right] = get_undirected_graph(conn_table, selection_cond, ...
    source_data_details, connectome_cond)

% Gets the IDs of presynaptic and postsynaptic neurons
pre_synaptic = unique([table2array(unique(conn_table(:,1))); ...
    table2array(unique(conn_table(:,2)))]);
post_synaptic = pre_synaptic;

% Computes the adjacency matrix
[adj_matrix, curr_name, curr_ID] = get_adjacency_matrix(conn_table, pre_synaptic, post_synaptic);

% Sorts the adjacency matrix according to hemisphere, layer and celltype
[adj_matrix_sorted, curr_name, idx_left, idx_right, a_left, a_right] = ...
    get_sorted_adjaceny_matrix(adj_matrix, curr_name, source_data_details, connectome_cond);

% Computes a symmetrix adjacency matrix
AT(:,:,1) = adj_matrix;
AT(:,:,2) = adj_matrix';
A = mean(AT,3);

% Computes the edges of an undirected graph (entries of symmetric adjacency matrix)
% and, in same loop, sets autapses of directed graph to NaN
undirected_graph_edges = [];
curr_adj_matrix = adj_matrix;
for loop_idx = 1:size(A,1)
    undirected_graph_edges = [undirected_graph_edges, A(loop_idx,loop_idx+1:end)];
    curr_adj_matrix(loop_idx,loop_idx) = NaN;
end
undirected_graph_edges = undirected_graph_edges';

% Computes the edges of a directed graph (entries of original adjacency matrix)
directed_graph_edges = reshape(curr_adj_matrix,size(curr_adj_matrix,1)*size(curr_adj_matrix,2),1);
directed_graph_edges = directed_graph_edges(~isnan(directed_graph_edges(:,1)),1);

end

function [adj_matrix, curr_name, curr_ID] = get_adjacency_matrix(curr_table, pre_synaptic, post_synaptic)

% Computes adjacency matrix

curr_name = [];
for loop_idx1 = 1:length(pre_synaptic)

    % Takes all connections in table that beliong to current presynaptic neuron
    loop_array1 = table2array(curr_table(find(table2array(curr_table(:,1)) == ...
        pre_synaptic(loop_idx1)),1:3));
    
    % Keeps track of presynaptic names
    curr_name = [curr_name; unique(curr_table(find(table2array(curr_table(:,1)) == ...
        pre_synaptic(loop_idx1)),4))];

    for loop_idx2 = 1:length(post_synaptic)
        
        % For the current presynaptic neuron, loops through postsynaptic
        % partners
        loop_array2 = loop_array1((loop_array1(:,2)) == post_synaptic(loop_idx2),:);
            
        if isempty(loop_array2) == 1
            % If there is no connection, a zero is placed into the adjacency matrix
            adj_matrix(loop_idx1,loop_idx2) = 0;
        else

            % If there is (at least) one connection, the n of synapses is placed into the adjacency matrix
            % Takes the sum over connections between two neurons: in MB plus CRE, for example)
            % If neurons connect in only one neuropil: it just takes the number of synapses there
            if size(loop_array2,1) > 1
                if isstring(loop_array2) == 1
                    % Optional: warning('!More than one row! If unexpected, check connectivity table!')
                    adj_matrix(loop_idx1,loop_idx2) = sum(str2double(loop_array2(:,3)));

                else
                    % Optional: warning('!More than one row! If unexpected, check connectivity table!')
                    adj_matrix(loop_idx1,loop_idx2) = sum(loop_array2(:,3));
                end
                
            else
                if isstring(loop_array2) == 1
                    adj_matrix(loop_idx1,loop_idx2) = str2double(loop_array2(1,3));
                else
                    adj_matrix(loop_idx1,loop_idx2) = loop_array2(1,3);
                end
            end
                     
        end

    end
end

curr_ID = pre_synaptic;

end

function symmetry_and_density_measures = get_symmetry_and_density_measures(adj_matrix, source_data_details, connectome_cond)

% Sets autapses in adjancency matrix to NaN
NaN_diag_AM = adj_matrix;
for loop_idx1 = 1:size(adj_matrix,1)
    NaN_diag_AM(loop_idx1,loop_idx1) = NaN;
end

% Computes correlation coefficient between adjacency matrix and its transpose
NaN_AM_transpose = NaN_diag_AM.';
symmetry_and_density_measures.correlation_with_transpose = ...
    corrcoef(NaN_AM_transpose, NaN_diag_AM, 'Rows', 'complete');
RS_NaN_diag_AM = reshape(NaN_diag_AM,...
    size(NaN_diag_AM,1)*size(NaN_diag_AM,2),1);

% Computes density without autapses
RS_NaN_diag_AM = RS_NaN_diag_AM(~isnan(RS_NaN_diag_AM),1);
symmetry_and_density_measures.density_adj = nnz(RS_NaN_diag_AM)/numel(RS_NaN_diag_AM);

% Quantifies reciprocity
adj_matrix_1 = adj_matrix;
threshold_value = 0:1:500;
reciprocity_measure = [];
for loop_idx2 = 1:length(threshold_value)
    curr_adj_matrix = double(adj_matrix_1 > threshold_value(loop_idx2));
    reciprocity_measure(loop_idx2) = (nnz(curr_adj_matrix.*(curr_adj_matrix.')))/((nnz(curr_adj_matrix)));
    n_connections(loop_idx2) = nnz(curr_adj_matrix);
end

symmetry_and_density_measures.reciprocity_measure = reciprocity_measure;
symmetry_and_density_measures.n_connections = n_connections;
symmetry_and_density_measures.threshold_value = threshold_value;

% Numerical estimate of symmetry of adjacency matrix using the squared
% Frobenius norm (i.e., the sum of squares)
% First, generates nullmodels by permuting the adjacency matrix and
% computes asymmetry of those nullmodels.
num_nullmodels = 100000;
norm_null_1 = NaN(num_nullmodels,1);
for loop_idx3 = 1:num_nullmodels
    curr_matrix = reshape(adj_matrix_1,size(adj_matrix_1,1)*size(adj_matrix_1,2),1);
    curr_matrix = curr_matrix(randperm(length(curr_matrix)),1);
    curr_matrix = reshape(curr_matrix,size(adj_matrix_1,1),size(adj_matrix_1,2));
    norm_null_1(loop_idx3,1) = (norm(((curr_matrix.')-curr_matrix),"fro").^2)/(norm(curr_matrix,"fro").^2); 
end

symmetry_and_density_measures.norm_null_1 = norm_null_1;

% Then computes asymmetry of actual adjacency matrix.
symmetry_and_density_measures.norm_actual = ...
    (norm(((adj_matrix_1.')-adj_matrix_1),"fro").^2)/(norm(adj_matrix_1,"fro").^2); 

if isempty(source_data_details) == 0

     if connectome_cond == "hemibrain"
         curr_file_name = source_data_details.file_name_hemibrain;
     elseif connectome_cond == "flywire"
         curr_file_name = source_data_details.file_name_flywire;    
     elseif connectome_cond == "MB"
         curr_file_name = source_data_details.file_name_MB;
     end
    
    reciprocity_table = table;
    reciprocity_table.reciprocity_measure = reciprocity_measure.';
    reciprocity_table.n_connections = n_connections.';
    reciprocity_table.threshold_value = threshold_value.';
    writetable(reciprocity_table, [source_data_details.data_path, curr_file_name],...
        'Sheet', 'reciprocity')

    asymmetry_table = table;
    asymmetry_table.asymmetry_connectivity_matrix = NaN(num_nullmodels,1);
    asymmetry_table.asymmetry_connectivity_matrix(1) = symmetry_and_density_measures.norm_actual;
    asymmetry_table.asymmetry_null_models = norm_null_1;
    writetable(asymmetry_table, [source_data_details.data_path, curr_file_name],...
        'Sheet', 'asymmetry')

end


end

function [ipsi_weights, contra_weights, size_ipsi, size_contra] = get_conn_weights(conn_table, dFB_layer_pre,...
    dFB_layer_post, cell_type_cond, consider_sub_layers)


% First, considers connections starting from left hemisphere.
if consider_sub_layers == "yes"
    curr_conn_table_left_to_right = conn_table(contains(table2array(conn_table(:,4)),'L') &...
        contains(table2array(conn_table(:,4)),string(dFB_layer_pre)) & ...
        contains(table2array(conn_table(:,5)),'R') & contains(table2array(conn_table(:,5)),string(dFB_layer_post)),:);
    curr_conn_table_left_to_left = conn_table(contains(table2array(conn_table(:,4)),'L') & ...
        contains(table2array(conn_table(:,4)),string(dFB_layer_pre)) & ...
        contains(table2array(conn_table(:,5)),'L') & contains(table2array(conn_table(:,5)),string(dFB_layer_post)),:);
else
    curr_conn_table_left_to_right = conn_table(contains(table2array(conn_table(:,4)),'L')  & ...
        contains(table2array(conn_table(:,5)),'R'),:);
    curr_conn_table_left_to_left = conn_table(contains(table2array(conn_table(:,4)),'L') & ...
        contains(table2array(conn_table(:,5)),'L'),:);
end


if cell_type_cond == "within" 
    conn_table_left_to_right = curr_conn_table_left_to_right(string(table2array(curr_conn_table_left_to_right(:,6))) == ...
        string(table2array(curr_conn_table_left_to_right(:,7))),:);
    conn_table_left_to_left = curr_conn_table_left_to_left(string(table2array(curr_conn_table_left_to_left(:,6))) == ...
        string(table2array(curr_conn_table_left_to_left(:,7))),:);
elseif cell_type_cond == "between"
    conn_table_left_to_right = curr_conn_table_left_to_right(string(table2array(curr_conn_table_left_to_right(:,6))) ~= ...
        string(table2array(curr_conn_table_left_to_right(:,7))),:);
    conn_table_left_to_left = curr_conn_table_left_to_left(string(table2array(curr_conn_table_left_to_left(:,6)))  ~= ...
        string(table2array(curr_conn_table_left_to_left(:,7))),:);
elseif cell_type_cond == "all"
    conn_table_left_to_right = curr_conn_table_left_to_right;
    conn_table_left_to_left = curr_conn_table_left_to_left;
else
    disp("!!!WARNING!!! WRONG CELL_TYPE_COND")
end


% Then, considers connections starting from right hemisphere.
if consider_sub_layers == "yes"
    curr_conn_table_right_to_left = conn_table(contains(table2array(conn_table(:,4)),'R') & ...
        contains(table2array(conn_table(:,4)),string(dFB_layer_pre)) & ...
        contains(table2array(conn_table(:,5)),'L') & contains(table2array(conn_table(:,5)),string(dFB_layer_post)),:);
    curr_conn_table_right_to_right = conn_table(contains(table2array(conn_table(:,4)),'R') & ...
        contains(table2array(conn_table(:,4)),string(dFB_layer_pre)) & ...
        contains(table2array(conn_table(:,5)),'R') & contains(table2array(conn_table(:,5)),string(dFB_layer_post)),:);
else       
    curr_conn_table_right_to_left = conn_table(contains(table2array(conn_table(:,4)),'R')  & ...
    contains(table2array(conn_table(:,5)),'L'),:);
    curr_conn_table_right_to_right = conn_table(contains(table2array(conn_table(:,4)),'R') & ...
    contains(table2array(conn_table(:,5)),'R'),:); 
end

if cell_type_cond == "within" 
    conn_table_right_to_left = curr_conn_table_right_to_left(string(table2array(curr_conn_table_right_to_left(:,6))) == ...
        string(table2array(curr_conn_table_right_to_left(:,7))),:);
    conn_table_right_to_right = curr_conn_table_right_to_right(string(table2array(curr_conn_table_right_to_right(:,6))) == ...
        string(table2array(curr_conn_table_right_to_right(:,7))),:);
elseif cell_type_cond == "between"
    conn_table_right_to_left = curr_conn_table_right_to_left(string(table2array(curr_conn_table_right_to_left(:,6))) ~= ...
        string(table2array(curr_conn_table_right_to_left(:,7))),:);
    conn_table_right_to_right = curr_conn_table_right_to_right(string(table2array(curr_conn_table_right_to_right(:,6))) ~= ...
        string(table2array(curr_conn_table_right_to_right(:,7))),:);
elseif cell_type_cond == "all"
    conn_table_right_to_left = curr_conn_table_right_to_left;
    conn_table_right_to_right = curr_conn_table_right_to_right;
else
    disp("!!!WARNING!!! WRONG CELL_TYPE_COND")
end

ipsi_weights = [table2array(conn_table_left_to_left(:,3)); table2array(conn_table_right_to_right(:,3))];
size_ipsi = size(unique([table2array(conn_table_left_to_left(:,1)); table2array(conn_table_right_to_right(:,1))]));

contra_weights = [table2array(conn_table_left_to_right(:,3)); table2array(conn_table_right_to_left(:,3))];
size_contra = size(unique([table2array(conn_table_left_to_right(:,1)); table2array(conn_table_right_to_left(:,1))]));

end

function [adj_matrix_sorted, curr_name, idx_left, idx_right, a_left, a_right] = get_sorted_adjaceny_matrix(...
    adj_matrix, curr_name, source_data_details, connectome_cond)

% Gets cell names and their original ordering
name_character_array = char(curr_name{:,:});
orig_idx = 1:size(curr_name,1);

% Finds cell names of left hemisphere
left_logical = contains(string(name_character_array),'_L');
left_names = curr_name(left_logical,:);
idx_left = orig_idx(left_logical);
[a_left, b_left] = sortrows(left_names);
idx_left = idx_left(b_left);

% Finds cell names of right hemisphere
right_logical = contains(string(name_character_array),'_R');
right_names = curr_name(right_logical,:);
idx_right = orig_idx(right_logical);
[a_right, b_right] = sortrows(right_names);
idx_right = idx_right(b_right);

% Sorts adjacency matrix
b = [idx_left, idx_right];
adj_matrix_sorted = adj_matrix(b, b);
curr_name = curr_name(b,1);


if isempty(source_data_details) == 0

     if connectome_cond == "hemibrain"
         curr_file_name = source_data_details.file_name_hemibrain;
     elseif connectome_cond == "flywire"
         curr_file_name = source_data_details.file_name_flywire;    
     elseif connectome_cond == "MB"
         curr_file_name = source_data_details.file_name_MB;
     end
    
    connectivity_matrix_cell_names = ([a_left; a_right]);
    writetable(connectivity_matrix_cell_names, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'cell_names')
    connectivity_matrix = table;
    connectivity_matrix.connectivity_matrix_sorted = adj_matrix_sorted;
    writetable(connectivity_matrix, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'connectivity_matrix')

end


end

function [conn_table_flywire, details_array_nu] = get_FlyWire_conn_table(source_data_details)


load('full_connections_table.mat');
classification_table = readtable("classification.csv");
load('classification_precise_ids.mat');
nu_class = table;
full_thing = table;
for loop_idx1 = 1:size(classification_precise_IDs,1)
    yo = classification_precise_IDs(loop_idx1,:);
    nu_class.root_id = classification_precise_IDs(loop_idx1,:);
    full_thing = [full_thing; nu_class];
end
classification_table = [full_thing, classification_table(:,2:end)];



accession_codes = find(contains(string(classification_table.cell_type),'FB6E1') | ...
    contains(string(classification_table.cell_type),'FB6E3') | ...
    contains(string(classification_table.cell_type),'FB6E4') | ...
    contains(string(classification_table.cell_type),'FB7E2'));
curr_names = classification_table(accession_codes,:);
full_IDs = classification_precise_IDs(accession_codes,:);

id_pre = full_connections_table.pre_pt_root_id;
id_post = full_connections_table.post_pt_root_id;
member_ids = ismember(id_pre, table2array(curr_names(:,1))) & ...
    ismember(id_post, table2array(curr_names(:,1)));
conn_table = full_connections_table(member_ids,:);

restrict_to_FB = 0;
if restrict_to_FB == 1 % condition if dFBN-to-dFBN connections outside the FB are considered
    conn_table = conn_table(contains(string(conn_table.neuropil), "FB"),:);
end
curr_hemi_IDs = classification_table(ismember(classification_table.root_id, table2array(curr_names(:,1))),:); 
curr_hemi_IDs = curr_hemi_IDs(curr_hemi_IDs.super_class == "central",:);
details_array = curr_hemi_IDs(:,["root_id","cell_type", "hemibrain_type","side"]);
details_array(contains(string(details_array.cell_type),'FB') == 0 & ...
    contains(string(details_array.hemibrain_type),'FB') == 0,:) = [];

details_array_nu = table;
details_array_nu.full_IDs = full_IDs((full_IDs) == details_array.root_id); % new
details_array_nu.cell_type = details_array.cell_type;
details_array_nu.hemibrain_type = details_array.hemibrain_type;



for loop_idx2 = 1:size(details_array,1)

    if details_array.side(loop_idx2,1) == "left"
        string_1 = string(details_array.cell_type(loop_idx2,:));
        string_2 = "_L";
        curr_sided_name = append(string_1, string_2);
        sided_name_array(loop_idx2,1) = curr_sided_name;
    elseif details_array.side(loop_idx2,1) == "right"
        string_1 = string(details_array.cell_type(loop_idx2,:));
        string_2 = "_R";
        curr_sided_name = append(string_1, string_2);
        sided_name_array(loop_idx2,1) = curr_sided_name;
    end

end

details_array_sided = details_array;
details_array_sided.cell_type = sided_name_array;
conn_table_flywire = conn_table(:,[1,2,4]);

for loop_idx3 = 1:size(conn_table_flywire,1)
    idx_array_start = table2array(details_array_sided(:,1)) == table2array(conn_table_flywire(loop_idx3,1));
    start_ID_sided(loop_idx3,1) = details_array_sided.cell_type(idx_array_start,1);
    start_ID_general(loop_idx3,1) = details_array.cell_type(idx_array_start,1);
    idx_array_end = table2array(details_array_sided(:,1)) == table2array(conn_table_flywire(loop_idx3,2));
    end_ID_sided(loop_idx3,1) = details_array_sided.cell_type(idx_array_end,1);
    end_ID_general(loop_idx3,1) = details_array.cell_type(idx_array_end,1);
end
 
conn_table_flywire = [conn_table_flywire, array2table(start_ID_sided), ...
    array2table(end_ID_sided), array2table(start_ID_general), ...
    array2table(end_ID_general)];




if isempty(source_data_details) == 0
    

    curr_file_name = source_data_details.file_name_flywire;
    accession_codes_and_names = details_array_nu;
    writetable(accession_codes_and_names, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'accession_codes_and_names')
    
    individual_connections = conn_table_flywire;
    writetable(individual_connections, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'individual_connections')
    
end


end

function [conn_table_MB_flywire, adj_matrix_MB, symmetry_and_density_measures_MB] = get_MB_conn_table(source_data_details)


load('full_connections_table.mat');
classification_table = readtable("classification.csv");
load('classification_precise_ids.mat');
nu_class = table;
full_t = table;
for loop_idx1 = 1:size(classification_precise_IDs,1)
    nu_class.root_id = classification_precise_IDs(loop_idx1,:);
    full_t = [full_t; nu_class];
end
classification_table = [full_t, classification_table(:,2:end)];


cell_names = "MBON";
name_cond = "MBON";
accession_codes = find(contains(string(classification_table.cell_type),cell_names) |...
    contains(string(classification_table.hemibrain_type),cell_names));
curr_names = classification_table(accession_codes,:);
full_IDs = classification_precise_IDs(accession_codes,:);

id_pre = full_connections_table.pre_pt_root_id;
id_post = full_connections_table.post_pt_root_id;
member_ids = ismember((id_pre), table2array(curr_names(:,1))) & ...
    ismember((id_post), table2array(curr_names(:,1)));

conn_table = full_connections_table(member_ids,:);

curr_hemi_IDs = classification_table(ismember((classification_table.root_id), ...
    table2array(curr_names(:,1))),:);
curr_hemi_IDs = curr_hemi_IDs(curr_hemi_IDs.super_class == "central",:);
details_array = curr_hemi_IDs(:,["root_id","cell_type", "hemibrain_type","side"]);
details_array(contains(string(details_array.cell_type),string(name_cond)) == 0 & ...
    contains(string(details_array.hemibrain_type),string(name_cond)) == 0,:) = [];

% Checks if flywire has assigned a name to cell, otherwise use hemibrain type
for loop_idx2 = 1:size(details_array,1)
    if size(char(details_array.cell_type(loop_idx2,1)),2) == 0
        details_array.cell_type(loop_idx2) = details_array.hemibrain_type(loop_idx2);
    end
end

details_array_nu = table;
details_array_nu.full_IDs = full_IDs((full_IDs) == details_array.root_id);
details_array_nu.cell_type = details_array.cell_type;
details_array_nu.hemibrain_type = details_array.hemibrain_type;

details_array_sided = details_array;
for loop_idx3 = 1:size(details_array,1)

    if details_array.side(loop_idx3,1) == "left"
        string_1 = string(details_array.cell_type(loop_idx3,:));
        string_2 = "_L";
        name_sided = append(string_1, string_2);
        sided_name_array(loop_idx3,1) = name_sided;
    elseif details_array.side(loop_idx3,1) == "right"
        string_1 = string(details_array.cell_type(loop_idx3,:));
        string_2 = "_R";
        name_sided = append(string_1, string_2);
        sided_name_array(loop_idx3,1) = name_sided;
    end

end

details_array_sided.cell_type = sided_name_array;

conn_table_MB_flywire = conn_table(:,[1,2,4]);

for loop_idx4 = 1:size(conn_table_MB_flywire,1)
    idx_array_start = table2array(details_array_sided(:,1)) == table2array(conn_table_MB_flywire(loop_idx4,1));
    start_ID_sided(loop_idx4,1) = details_array_sided.cell_type(idx_array_start,1);
    start_ID_general(loop_idx4,1) = details_array.cell_type(idx_array_start,1);
    idx_array_end = table2array(details_array_sided(:,1)) == table2array(conn_table_MB_flywire(loop_idx4,2));
    end_ID_sided(loop_idx4,1) = details_array_sided.cell_type(idx_array_end,1);
    end_ID_general(loop_idx4,1) = details_array.cell_type(idx_array_end,1);
end
 
conn_table_MB_flywire = [conn_table_MB_flywire, array2table(start_ID_sided), ...
    array2table(end_ID_sided), array2table(start_ID_general), ...
    array2table(end_ID_general)];

[directed_table, undirected_table, ~, adj_matrix_MB, adj_matrix_sorted, ~, ~, ~, ~, ~, ~] = ...
    get_undirected_graph(conn_table_MB_flywire, 'all', ...
    source_data_details, "MB");
symmetry_and_density_measures_MB = get_symmetry_and_density_measures(adj_matrix_MB, source_data_details, "MB");
symmetry_and_density_measures_MB.adj_matrix_sorted = adj_matrix_sorted;
symmetry_and_density_measures_MB.details_array_nu = details_array_nu;



get_fig_panel_log_log_dist_v2(directed_table, undirected_table, get_color, 3, 3, source_data_details, "MB")





if isempty(source_data_details) == 0
    


    curr_file_name = source_data_details.file_name_MB;
    accession_codes_and_names = details_array_nu;
    writetable(accession_codes_and_names, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'accession_codes_and_names')

    individual_connections = conn_table_MB_flywire;
    writetable(individual_connections, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'individual_connections')
    
end



end

%% Plotting functions

function get_fig_panel_log_log_dist_v2(directed_table, undirected_table, color, x_pos, y_pos,...
    source_data_details, connectome_cond)


if connectome_cond == "hemibrain"
    n_bins = 12;
elseif connectome_cond == "flywire"
    n_bins = 11;
elseif connectome_cond == "MB"
    n_bins = 7;
end

curr_table_directed = log10(directed_table);
curr_table_undirected = log10(undirected_table);

cond_include_zeros = 0;
if cond_include_zeros == 1
    curr_table_directed(curr_table_directed == -inf,1) = 0;
    curr_table_undirected(curr_table_undirected == -inf,1) = 0;
elseif cond_include_zeros == 0
    curr_table_directed = curr_table_directed(curr_table_directed ~= -inf,1);
    curr_table_undirected = curr_table_undirected(curr_table_undirected ~= -inf,1);
end

[n_counts_directed, n_synapses_directed] = ...
    histcounts(curr_table_directed, 'NumBins', n_bins, 'BinLimits', [0, 3],...
    'Normalization', 'probability');
[n_counts_undirected, n_synapses_undirected] = ...
    histcounts(curr_table_undirected, 'NumBins', n_bins, 'BinLimits', [0, 3],...
    'Normalization', 'probability');

n_synapses_directed = n_synapses_directed-(mean(diff(n_synapses_directed))/2);
n_synapses_directed = n_synapses_directed(2:end);
n_synapses_undirected = n_synapses_undirected-(mean(diff(n_synapses_undirected))/2);
n_synapses_undirected = n_synapses_undirected(2:end);

sz_1 = 1.8;
ht_1 = 1.8;
xlm_1 = [0, 3];
ylm_1 = [-3, 0];
major_x_ticks = 1;
major_y_ticks = 1.5;
plotting_color_directed = color.dark_gray;
plotting_color_undirected = color.light_gray;

pnl_1 = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos, sz_1, ht_1]);
hold on
plot(pnl_1, n_synapses_undirected, log10(n_counts_undirected),...
    'Color', plotting_color_undirected, 'LineWidth', 1.5)
scatter(pnl_1, n_synapses_directed, log10(n_counts_directed), 10,...
    'MarkerFaceColor', plotting_color_directed, 'MarkerEdgeColor', 'none')
get_default_separated_ax(pnl_1, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ...
    ylm_1(1), ylm_1(2), ylm_1(1), ylm_1(2), major_x_ticks, major_y_ticks,...
    "linear", "linear", 'Log_{10}({\itn} synapses)',...
    {'Log_{10}';'(fraction of connections)'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)


if isempty(source_data_details) == 0

     if connectome_cond == "hemibrain"
         curr_file_name = source_data_details.file_name_hemibrain;
     elseif connectome_cond == "flywire"
         curr_file_name = source_data_details.file_name_flywire;    
     elseif connectome_cond == "MB"
         curr_file_name = source_data_details.file_name_MB;
     end
    
    directed_vs_undirected = table;
    directed_vs_undirected.n_synapses_undirected = n_synapses_undirected.';
    directed_vs_undirected.n_counts_undirected = n_counts_undirected.';
    directed_vs_undirected.n_synapses_directed = n_synapses_directed.';
    directed_vs_undirected.n_counts_directed = n_counts_directed.';
    writetable(directed_vs_undirected, ...
        [source_data_details.data_path, curr_file_name],...
        'Sheet', 'directed_vs_undirected')

end


end

function get_fig_panel_weight_histogram(curr_weights_1, curr_weights_2, ...
    x_pos, y_pos, color, ylm_1, major_y_ticks, text_ipsi, text_contra,...
    text_headline, plot_name, source_data_details, compare_cond, connectome_cond)

if connectome_cond == "hemibrain"
    bin_lims = [0, 250];
    major_x_ticks = 125;
    xlm_1 = bin_lims;
    if isempty(source_data_details) == 0
        curr_filename = source_data_details.file_name_hemibrain;
    end
elseif connectome_cond == "flywire"
    bin_lims = [0, 150];
    major_x_ticks = 75;
    xlm_1 = bin_lims;
    if isempty(source_data_details) == 0
        curr_filename = source_data_details.file_name_flywire;
    end
end
n_bins = 30;
sz_1 = 1.8;
ht_1 = 1.8;

if compare_cond == "hemispheres"
    plotting_color_1 = color.light_gray;
    plotting_color_2 = color.navy;
    n_bins = 30;
    alpha_1 = 1;
    alpha_2 = 0.5;
    sheet_names = {'Within_hemisphere','Inter_hemisphere'};
elseif compare_cond == "layers"
    plotting_color_1 = color.light_gray;
    plotting_color_2 = [0, 0, 0];
    n_bins = 50;
    alpha_1 = 1;
    alpha_2 = 0.6;
    sheet_names = {'Between_layers','Within_layers'};
end

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
histogram(pnl_1, curr_weights_1, 'Normalization', 'probability',...
    'NumBins', n_bins, 'BinLimits', bin_lims, 'EdgeColor', 'none', 'FaceColor',...
    plotting_color_1, 'FaceAlpha', alpha_1)    
hold on    
histogram(pnl_1, curr_weights_2, 'Normalization', 'probability',...
    'NumBins', n_bins, 'BinLimits', bin_lims, 'EdgeColor', 'none', 'FaceColor',...
    plotting_color_2, 'FaceAlpha', alpha_2)
get_default_separated_ax(pnl_1, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ...
    ylm_1(1), ylm_1(2), ylm_1(1), ylm_1(2),...
    major_x_ticks, major_y_ticks, "linear", "linear",...
    '{\itn} synapses', string(text_headline),...
    'k', 'k', 'k', 'k', sz_1, ht_1)

text_pos = x_pos+sz_1;
y_pos_ipsi = y_pos+ht_1+0.2;
if compare_cond == "hemispheres"
    get_default_annotation(text_pos, y_pos_ipsi, 'Within hemispheres', plotting_color_1, 'normal', "right")
    get_default_annotation(text_pos, y_pos_ipsi+0.25, 'Between hemispheres', plotting_color_2, 'normal', "right")
elseif compare_cond == "layers"
    get_default_annotation(text_pos, y_pos_ipsi, 'Between layers', plotting_color_1, 'normal', "right")
    get_default_annotation(text_pos, y_pos_ipsi+0.25, 'Within same layer', plotting_color_2, 'normal', "right")
end

if isempty(source_data_details) == 0
    % Saves the source data
    n_syn = NaN(1000,2);
    n_syn(1:length(curr_weights_1),1) = curr_weights_1';
    n_syn(1:length(curr_weights_2),2) = curr_weights_2';
    n_syn = array2table(n_syn);
    n_syn.Properties.VariableNames = string(sheet_names);
    writetable(n_syn, ...
        [source_data_details.data_path, curr_filename],'Sheet',string(plot_name))
end

end

function get_fig_panel_adjacency_matrix_heatmap(x_pos, y_pos, adj_matrix, idx_left, idx_right, ...
    a_left, a_right, color, connectome_cond)


s_z = 5;
xlm_1 = [0.5, size(adj_matrix,1)+0.5];
ylm_1 = [0.5, size(adj_matrix,1)+0.5];
clim_adj = [0, max(max(adj_matrix))];
width_layer = 1.5;
width_hemi = 1.5;
color_layer6 = color.dark_yellow;
color_layer7 = color.medium_gray;

curr_a_left = char(a_left{:,:});
curr_a_left = curr_a_left(:,3:6);
n_6L = sum(contains(string(curr_a_left),"6"));

curr_a_right = char(a_right{:,:});
curr_a_right = curr_a_right(:,3:6);
n_6R = sum(contains(string(curr_a_right),"6"));
 
curr_colormap = get_default_color_map(5000, 1, [0, 0, 0], "normal");


%% Plots data
min_max = [min(adj_matrix,[],"all"), max(adj_matrix,[],"all")];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, s_z, s_z]);
imagesc(curr_pnl, 1:size(adj_matrix,2), 1:size(adj_matrix,2), adj_matrix') % transpose because rows are presypases
colormap(curr_pnl, curr_colormap)
clim(curr_pnl, clim_adj)
curr_pnl.YDir = "normal";
curr_pnl.Box = 'on';
curr_pnl.Color = 'none';
curr_pnl.TickLength = [0, 0];
curr_pnl.XLabel.Visible = 'off';
curr_pnl.YLabel.Visible = 'off';
curr_pnl.YLabel.Color = 'none';
curr_pnl.XTick = [];
curr_pnl.YTick = [];
get_default_colorbar(curr_pnl, x_pos+s_z+0.1, y_pos, 0.2, s_z,...
    clim_adj, min_max, '{\itn} synapses', 'eastoutside')

additional_spread = 0.4;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-0.2, s_z, 0.1]);
hold on
plot(curr_pnl, [1-additional_spread, length(idx_left)+additional_spread], [0, 0],'Color', 'k', 'LineWidth', width_hemi)
plot(curr_pnl, [1-additional_spread, n_6L+additional_spread], [1, 1], 'Color', color_layer6, 'LineWidth', width_layer)
plot(curr_pnl, [n_6L+1-additional_spread, length(idx_left)+additional_spread], [1, 1],...
    'Color', color_layer7, 'LineWidth', width_layer)

plot(curr_pnl, [length(idx_left)+1-additional_spread, size(adj_matrix,1)+additional_spread], [0, 0],...
    'Color', 'k', 'LineWidth', width_hemi)
plot(curr_pnl, [(length(idx_left)+1-additional_spread), (length(idx_left)+1)+(n_6R-1)+additional_spread], [1, 1],...
    'Color', color_layer6, 'LineWidth', width_layer)
plot(curr_pnl, [(length(idx_left)+1-additional_spread)+(n_6R-1)+1, size(adj_matrix,1)+additional_spread], [1, 1],...
    'Color', color_layer7, 'LineWidth', width_layer)
get_default_ax(curr_pnl, xlm_1(1), xlm_1(2), [], [], 0, 1,...
    [], [], [], [],"linear", "linear", [], [], 'none', 'none', 'none', 'none', s_z, 0.1)


curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos-0.2, y_pos, 0.1, s_z]);
hold on
plot(curr_pnl, [0, 0], [1-additional_spread, length(idx_left)+additional_spread], 'Color', 'k', 'LineWidth', width_hemi)
plot(curr_pnl, [1, 1], [1-additional_spread, n_6L+additional_spread], 'Color', color_layer6, 'LineWidth', width_layer)
plot(curr_pnl, [1, 1], [n_6L+1-additional_spread, length(idx_left)+additional_spread],...
    'Color', color_layer7, 'LineWidth', width_layer)

plot(curr_pnl, [0, 0], [length(idx_left)+1-additional_spread, size(adj_matrix,1)+additional_spread],...
    'Color', 'k', 'LineWidth', width_hemi)
plot(curr_pnl, [1, 1], [(length(idx_left)+1-additional_spread), (length(idx_left)+1)+(n_6R-1)+additional_spread],...
    'Color', color_layer6, 'LineWidth', width_layer)
plot(curr_pnl, [1, 1], [(length(idx_left)+1-additional_spread)+(n_6R-1)+1, size(adj_matrix,1)+additional_spread],...
    'Color', color_layer7, 'LineWidth', width_layer)
get_default_ax(curr_pnl, 0, 1, [], [], ylm_1(1), ylm_1(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.1, s_z)

get_default_annotation(x_pos+(s_z/4), y_pos-0.5, 'Left hemisphere', 'k', 'normal', "center")
get_default_annotation(x_pos+((s_z/4)*3), y_pos-0.5, 'Right hemisphere', 'k', 'normal', "center")
get_default_annotation_rotated(x_pos-0.5, y_pos+(s_z/4)-0.5, 'Left hemisphere', 'k', 'normal', "left")
get_default_annotation_rotated(x_pos-0.5, y_pos+((s_z/4)*3)-0.5, 'Right hemisphere', 'k', 'normal', "left")
get_default_annotation(x_pos+s_z-1.1, y_pos+s_z+0.3, 'Layer 6', color_layer6, 'normal', "right")
get_default_annotation(x_pos+s_z, y_pos+s_z+0.3, 'Layer 7', color_layer7, 'normal', "right")

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos-0.5, y_pos-0.5, 0.5, 0.5]);
hold on
plot(curr_pnl, [0, 1], [0, 0], 'Color', 'k', 'LineWidth', get_default_scale_bar_width)
plot(curr_pnl, [0, 0], [0, 1], 'Color', 'k', 'LineWidth', get_default_scale_bar_width)
get_default_ax(curr_pnl, 0, 1, [], [], 0, 1, [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.1, s_z)
get_default_annotation(x_pos-0.45, y_pos-0.7, 'Pre', 'k', 'normal', "left")
get_default_annotation_rotated(x_pos-0.8, y_pos-0.25, 'Post', 'k', 'normal', "left")


y_pos_annotation = y_pos-0.25;
if connectome_cond == "hemibrain"

    annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos, y_pos_annotation+s_z+0.3, 10, 0.25], ...
        'string', 'Hemibrain', 'EdgeColor', 'none', 'FontSize', 6, 'FontWeight', 'bold', 'Color', [0, 0, 0],...
        'FontAngle', 'normal', 'Margin', 0, 'HorizontalAlignment', "left")

elseif connectome_cond == "flywire"
    % FlyWire FAFB v783
    annotation('textbox', 'Units', 'centimeters', 'Position',[x_pos, y_pos_annotation+s_z+0.3, 10, 0.25], ...
        'string', 'FlyWire FAFB', 'EdgeColor', 'none', 'FontSize', 6, 'FontWeight', 'bold', 'Color', [0, 0, 0],...
        'FontAngle', 'normal', 'Margin', 0, 'HorizontalAlignment', "left")

end


end

function get_fig_panel_reciprocity(x_pos, y_pos, symmetry_and_density_measures_hemibrain,...
    symmetry_and_density_measures_flywire, symmetry_and_density_measures_MB, color)

sz_1 = 5.3;
ht_1 = 1.8;
ylm_1 = [0, 1];
xlm_1 = [0, 100];
major_x_ticks = 25;
major_y_ticks = 0.5;
color_mean = color.navy;
color_hemibrain = [0, 0, 0];
color_flywire = color.medium_gray;
color_MB = color.light_gray;

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
plot(pnl_1, symmetry_and_density_measures_hemibrain.threshold_value, ...
    symmetry_and_density_measures_hemibrain.reciprocity_measure,...
    'Color', color_hemibrain, "LineWidth", 1, 'LineStyle', '-')
hold on
plot(pnl_1, symmetry_and_density_measures_flywire.threshold_value, ...
    symmetry_and_density_measures_flywire.reciprocity_measure,...
    'Color', color_flywire, "LineWidth", 1, 'LineStyle', '-')
plot(pnl_1, symmetry_and_density_measures_hemibrain.threshold_value, ...
    mean([symmetry_and_density_measures_flywire.reciprocity_measure; ...
    symmetry_and_density_measures_hemibrain.reciprocity_measure],1),...
    'Color', color_mean, "LineWidth", 2, 'LineStyle', '-')
plot(pnl_1, symmetry_and_density_measures_MB.threshold_value, ...
    symmetry_and_density_measures_MB.reciprocity_measure,...
    'Color', color_MB, "LineWidth", 1, 'LineStyle', '-')
get_default_separated_ax(pnl_1, xlm_1(1), xlm_1(2), ...
    xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2), ylm_1(1), ylm_1(2),...
    major_x_ticks, major_y_ticks, "linear", "linear",...
    'Threshold ({\itn} synapses)', 'Reciprocity',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

text_pos = x_pos+sz_1;
y_dist_annotation = 0.2;
get_default_annotation(text_pos, y_pos+ht_1-y_dist_annotation, 'Hemibrain', color_hemibrain, 'normal', "right")
get_default_annotation(text_pos, y_pos+ht_1-y_dist_annotation-0.25, 'Flywire', color_flywire, 'normal', "right")
get_default_annotation(text_pos, y_pos+ht_1-y_dist_annotation-0.5, 'Mean', color_mean, 'normal', "right")
get_default_annotation(text_pos, y_pos+ht_1-y_dist_annotation-1.1, 'Mushroom Body (MBONs)', color_MB, 'normal', "right")

end

function get_fig_panel_num_symmetry(x_pos, y_pos, symmetry_and_density_measures_hemibrain,...
    symmetry_and_density_measures_flywire, symmetry_and_density_measures_MB, color)

sz_1 = 5.8;
ht_1 = 1.8;
ylm_1 = [0, 0.5];
xlm_1 = [0, 2.5];
major_x_ticks = 0.5;
major_y_ticks = 0.25;
color_hemibrain = [0, 0, 0];
color_flywire = color.medium_gray;
color_MB = color.light_gray;
color_symmetric = color.light_lavender;
height_stem = ylm_1(2)*0.75;
lw_hist = 0.5;
lw_stem = 0.5;
n_bins = 30;
m_size = 3;

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
histogram(pnl_1, symmetry_and_density_measures_MB.norm_null_1, ...
    'DisplayStyle', 'stairs', 'Normalization', 'probability', 'NumBins', n_bins, 'EdgeColor', color_MB,...
    "LineWidth", lw_hist)
histogram(pnl_1, symmetry_and_density_measures_hemibrain.norm_null_1, ...
    'DisplayStyle', 'stairs', 'Normalization', 'probability', 'NumBins', n_bins, 'EdgeColor', color_hemibrain,...
    "LineWidth", lw_hist)
histogram(pnl_1, symmetry_and_density_measures_flywire.norm_null_1, ...
    'DisplayStyle', 'stairs', 'Normalization', 'probability', 'NumBins', n_bins, 'EdgeColor', color_flywire,...
    "LineWidth", lw_hist)
stem(pnl_1, 0, height_stem,...
    'filled', 'MarkerSize', m_size, 'Color', color_symmetric, "LineWidth", lw_stem, 'LineStyle', '-')
stem(pnl_1, symmetry_and_density_measures_hemibrain.norm_actual, height_stem,...
    'filled', 'MarkerSize', m_size, 'Color', color_hemibrain, "LineWidth", lw_stem, 'LineStyle', '-')
stem(pnl_1, symmetry_and_density_measures_flywire.norm_actual, height_stem,...
    'filled', 'MarkerSize', m_size, 'Color', color_flywire, "LineWidth", lw_stem, 'LineStyle', '-')
stem(pnl_1, symmetry_and_density_measures_MB.norm_actual, height_stem,...
    'filled', 'MarkerSize', m_size, 'Color', color_MB, "LineWidth", lw_stem, 'LineStyle', '-')
get_default_separated_ax(pnl_1, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ...
    ylm_1(1), ylm_1(2), ylm_1(1), ylm_1(2), major_x_ticks, major_y_ticks, "linear", "linear",...
    'Asymmetry of connectivity matrix', {'Fraction';'(nullmodels)'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)

get_default_annotation(x_pos+sz_1, y_pos+ht_1, 'Hemibrain', color_hemibrain, 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos+ht_1-0.25, 'Flywire', color_flywire, 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos+ht_1-0.7, ["Symmetric";"matrix"], color_symmetric, 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos+0.25, 'MBONs', color_MB, 'normal', "right")

end
