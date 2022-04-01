function mutationChildren = mut_gaussiana_ind(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate)
global scale shrink lb ub
if(strcmpi(options.PopulationType,'doubleVector'))

    if nargin < 9 || isempty(shrink)
        shrink = 1;
    end
    if nargin < 8 || isempty(scale)
        scale = 1;
    end
    
    if (shrink > 1) || (shrink < 0)
        warning(message('globaloptim:mutationgaussian:shrinkFactor'));
    end

    scale = scale - shrink * scale * state.Generation/options.Generations;

    range = options.PopInitRange;
    lower = range(1,:);
    upper = range(2,:);
    scale = scale * (upper - lower);

    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        parent = thisPopulation(parents(i),:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r=randn(1,length(parent));
        rmax = max(r);
        rmin = min(r)-0.001;
        rpor = (r-rmin)*(1/(rmax-rmin));
        rf = rpor.*(ub-lb);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        mutationChildren(i,:) = parent  + scale .*rf;
    end
elseif(strcmpi(options.PopulationType,'bitString'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mutacao uniforme %%%%%%%%%%%%%%%%%%%%%%%    
if nargin < 8 || isempty(mutationRate)
    mutationRate = 0.01; % default mutation rate
end

if(strcmpi(options.PopulationType,'doubleVector'))
    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        child = thisPopulation(parents(i),:);
        % Each element of the genome has mutationRate chance of being mutated.
        mutationPoints = find(rand(1,length(child)) < mutationRate);
        % each gene is replaced with a value chosen randomly from the range.
        range = options.PopInitRange;
        % range can have one column or one for each gene.
        [r,c] = size(range);
        if(c ~= 1)
            range = range(:,mutationPoints);
        end   
        lower = range(1,:);
        upper = range(2,:);
        span = upper - lower;
        
        
        child(mutationPoints) = lower + rand(1,length(mutationPoints)) .* span;
        mutationChildren(i,:) = child;
    end
elseif(strcmpi(options.PopulationType,'bitString'))
    
    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        child = thisPopulation(parents(i),:);
        mutationPoints = find(rand(1,length(child)) < mutationRate);
        child(mutationPoints) = ~child(mutationPoints);
        mutationChildren(i,:) = child;
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
end