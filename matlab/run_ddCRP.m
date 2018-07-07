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

function run_ddCRP(surfacefile,labelfile,regions,datafile,outputfile,sizes,alpha,kappa,nu,sigsq,pass_limit)

    addpath('/mnt/parcellator/parcellation/Matlab/SurfaceDistanceComputations/');
    addpath('/mnt/parcellator/parcellation/Matlab/ddCRP/');
    
    fprintf('Loading surface file.\n');
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
        fprintf('Region: %s\n',char(r));
        indices = [indices;find(label.cdata == labelmap(char(r)))];
    end
    
    indices = sort(indices);
    
    fprintf('Filtering adjacency list.\n');
    filtered = filter_adjacency(adjacency_list,indices);
    
    fprintf('Loading data matrix.\n');
    try
        temp = load(datafile);
    catch
        Warning('Data file does not exist.');
    end
    fn = fieldnames(temp);
    data = temp.(fn{1});
    downsampled = data(indices,:);
    
    similarities = corr(downsampled');
    D = 1-similarities;
    D = normr(D);
    
    clear downsampled_data similarities adjacency_list
    
    [z,Z] = WardClustering(D,filtered,7);
    [map_z, stats] = InitializeAndRunddCRP(Z,D,filtered, sizes, alpha, kappa, nu, sigsq, pass_limit,[],true);
    
    outfunc = sprintf('%s.func.gii',outputfile);
    outstat = sprintf('%s.mat',outputfile);
    
    save(outfunc,'map_z','-v7.3');
    save(outstat,'stats','-v7.3');
    
end