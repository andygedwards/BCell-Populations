load_dir = '/Users/Andy/Backup/Drive - Sylloge Backup/Research/Pancreatic Islets/Manuscript/Code for publication/Glucose screen/results/Combined/3000/0.370/ka_0.00004';
low = load([load_dir,filesep,'low_glucose.mat']);
high = load([load_dir,filesep,'high_glucose.mat']);

pop_dir = '/Users/Andy/Backup/Drive - Sylloge Backup/Research/Pancreatic Islets/Manuscript/Code for publication/Create Population/Combined/3000/0.370';
load([pop_dir,filesep,'modParam_scaling_3000_0.370.mat'])
load([pop_dir,filesep,'modParam_names_3000_0.370.mat'])

low_mat = cell2mat(low.param_low);
high_mat = cell2mat(high.param_high);

[p,p_names] = Riz2014_init_parameters_INa_low;
modParam_scaling_new = modParam_scaling;
for i = 1:length(modParam_names)
    modParam_inds(i) = find(strcmp(p_names,modParam_names{i}));
    scaled_mat_low(:,i) = low_mat(:,modParam_inds(i))./p(modParam_inds(i));
    if strcmp(modParam_names{i},'g_Na')
       modParam_scaling_new(:,i) = scaled_mat_low(:,i);
    end
end
length(find(modParam_scaling_new-modParam_scaling > 0.00003))
modParam_scaling = modParam_scaling_new;
save([pop_dir,filesep,'modParam_scaling_3000_0.370.mat'],'modParam_scaling');






% Cellload_dir = '/Users/Andy/Backup/Drive - Sylloge Backup/Research/Pancreatic Islets/Manuscript/Code for publication/Glucose screen/data/Combined/3000/0.370/ka_0.0001/cells';
% % Cellload_dir_ref = '/Users/Andy/Backup/Drive - Sylloge Backup/Research/Pancreatic Islets/Manuscript/Code for publication/Glucose screen/data/Combined/3000/0.370/ka_0.0001 copy/cells';
% % nums = 1;
% nums = 1:1:3000;
% 
% [p,p_names] = Riz2014_init_parameters_INa_low;
% 
% count_low = 0;
% count_high = 0;
% for i = 1:length(nums)
%     cell_low = load([Cellload_dir,filesep,num2str(nums(i)),'_low.mat']);
%     cell_high = load([Cellload_dir,filesep,num2str(nums(i)),'_high.mat']);
%     % cell_low_ref = load([Cellload_dir_ref,filesep,num2str(nums(i)),'_low.mat']);
%     % cell_high_ref = load([Cellload_dir_ref,filesep,num2str(nums(i)),'_high.mat']);
%     % if ~isequal(cell_low.param,cell_low_ref.param)
%     %     inds_low(count_low+1) = i;
%     %     count_low = count_low+1;
%     % end
%     % if ~isequal(cell_high.param,cell_high_ref.param)
%     %     inds_high(count_high+1) = i;
%     %     count_high = count_high+1;
%     % end
%     low.param_low{i}(76) = 0.0001;
%     high.param_high{i}(76) = 0.0001;
%     if ~isequal(cell_low.param, low.param_low{i})
%         disp(['Low differences at ', num2str(i)])
%         find((low.param_low{i}-cell_low.param)~=0)
%         % cell_low.param = low.param_low{i};
%         % save_struct(cell_low,[Cellload_dir,filesep,num2str(nums(i)),'_low.mat'])
%         count_low = count_low+1;
%     end
%     cell_low.param = low.param_low{i};
%     save_struct(cell_low,[Cellload_dir,filesep,num2str(nums(i)),'_low.mat'])
%     if ~isequal(cell_high.param, high.param_high{i})
%         disp(['High differences at ', num2str(i)])
%         find((high.param_high{i}-cell_high.param)~=0);
%         % cell_high.param = high.param_high{i};
%         % save_struct(cell_high,[Cellload_dir,filesep,num2str(nums(i)),'_high.mat'])
%         count_high = count_high+1;
%     end
%     cell_high.param = high.param_high{i};
%     save_struct(cell_high,[Cellload_dir,filesep,num2str(nums(i)),'_high.mat'])
% end
% 
% count_low
% count_high
% 
% 
% 
% function save_struct(S,path)
%     fields = fieldnames(S);
%     % Assign each field to a variable
%     for i = 1:numel(fields)
%         eval([fields{i}, '= S.(fields{i});']);
%     end
%     % Save variables to a .mat file (overwriting any existing file)
%     save(path, fields{:});
% end
