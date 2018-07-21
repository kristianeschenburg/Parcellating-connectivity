% Wrapper function to run ddCRP method.
% Parameters:
% surfacefile: path to surface file
% labelfile: path to label file
% regions: ['region1','region2'] -- regions to include in adjacency matrix
% datafile: data from which to generate similarity matrix
% outputfile: basename of output data
% sizes: initialize expected number of clusters to generate
% alpha: hyperparameter (default = 10)
% kappa: hyperparameter (default = 0.0001)
% nu: hyperparameter (default = 1)
% sigsq: hyperparameter of cluster variance
% pass_limit: number of MCMC passes

function run_ddCRP(surfacefile, labelfile, datafile, outputfile, ... 
    sizes, alpha, kappa, nu, sigsq, pass_limit, ...
    varargin)

    addpath('./utilities/');

    % We add an optional parameter to receive the edge_prior matrix.
    p = inputParser;
    errorMsg = 'Value must be numeric or boolean.'; 
    validationFcn = @(x) assert(isnumeric(x) || islogical(x), errorMsg);
    p.addParameter('edge_prior',false,validationFcn);

    p.parse(varargin{:})
    edge_prior = p.Results.edge_prior;
    fprintf('Edge prior: \n');
    disp(edge_prior);

    if isa(sizes,'string') || isa(sizes,'char')
        sizes = str2double(sizes);
    end
    
    if isa(alpha,'string') || isa(alpha,'char')
        alpha = str2double(alpha);
    end
    
    if isa(kappa,'string') || isa(kappa,'char')
        kappa = str2double(kappa);
    end
    
    if isa(nu,'string') || isa(nu,'char')
        nu = str2double(nu);
    end
    
    if isa(sigsq,'string') || isa(sigsq,'char')
        sigsq = str2double(sigsq);
    end
    
    if isa(pass_limit,'string') || isa(pass_limit,'char')
        pass_limit = str2double(pass_limit);
    end

	fprintf('\nFile settings: \n');
    
    disp(surfacefile)
    
    [~,name,ext] = fileparts(surfacefile);
    testfile = sprintf('%s%s',name,ext);
	fprintf('Surface file: %s\n',testfile);
    
    [~,name,ext] = fileparts(labelfile);
    testfile = sprintf('%s%s',name,ext);
	fprintf('Label file: %s\n',testfile);
    
    [~,name,ext] = fileparts(datafile);
    testfile = sprintf('%s%s',name,ext);
	fprintf('Data file: %s\n',testfile);
    
    [~,name,ext] = fileparts(outputfile);
    testfile = sprintf('%s%s',name,ext);
	fprintf('Output File: %s\n',testfile);

	fprintf('\n');
   
	fprintf('Parameter settings: \n');
    fprintf('Sizes: %.3f\n',sizes);
	fprintf('Alpha: %.3f\n',alpha);
	fprintf('Kappa: %.3f\n',kappa);
	fprintf('Nu: %.3f\n',nu);
	fprintf('SigSQ: %.3f\n',sigsq);
	fprintf('Pass_Limit: %.3f\n',pass_limit);

    regions = {'L_inferiorparietal'; 'L_supramarginal'};

    fprintf('\nLoading surface file.\n');
    try
        surface = gifti(surfacefile);
    catch
        Warning('Surface file does not exist.');
    end
    
    fprintf('Generating full adjacency list.\n');
    adjacency_list = generate_adjacency(surface,false);
    
    fprintf('Loading label file.\n');
    try
        label = gifti(labelfile);
    catch
        Warning('Label file does not exist.');
    end
    
    labelmap = containers.Map(label.labels.name,label.labels.key);
    indices = [];
    for r = regions'
        indices = [indices;find(label.cdata == labelmap(char(r)))];
    end
    indices = sort(indices);
    
    
    fprintf('Filtering adjacency list.\n');
    filtered_adjacency = filter_adjacency(adjacency_list,indices);
    
    if edge_prior
        fprintf('Filtering edge prior matrix.\n');
        edge_prior = filter_prior(adjacency_list,indices,edge_prior);
    end
    
    fprintf('Loading data matrix.\n\n');
    try
        temp = load(datafile);
    catch
        Warning('Data file does not exist.');
    end
    
    fn = fieldnames(temp);
    data = temp.(fn{1});
    data = normr(data);
    downsampled = data(indices,:);
    
    S = corr(downsampled');
    S = S - diag(diag(S));
    
    S = normr(S);
    
    clear downsampled_data similarities adjacency_list
    
    [~,Z] = WardClustering(S,filtered_adjacency,sizes);
    [map_z, stats] = InitializeAndRunddCRP(Z, S, filtered_adjacency, ... 
        sizes, alpha, kappa, nu, sigsq, pass_limit, ... 
        [], true, 'edge_prior', edge_prior);
    
    outfunc = sprintf('%s.label.mat',outputfile);
    outstat = sprintf('%s.stats.mat',outputfile);
    
    map_z_full = zeros(size(label.cdata,1));
    map_z_full(indices) = map_z;
    
    save(outfunc,'map_z_full','-v7.3');
    save(outstat,'stats','-v7.3');
    
exit