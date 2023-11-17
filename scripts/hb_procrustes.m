function [E, E0, U_pt, S, S0] = hb_procrustes(U, varargin)
% HB_PROCRUSTES computes distance between two vector sets based Procrustes
% transform (PT). The vector sets are rotated to optimaly match each other.
% An ensemble measure of cosine similarity between non-matching vectors is
% returned, reflecting how different the two sets are.
%
% Inputs:
%   U: a cell array of length two. Each element is a matrix, where each
%   column represents an observation (e.g. an eigenmode). The size of the 
%   matrices should be the same.
%
% Outputs:
%   E: residual error after PT. 
%   E0: initial residual error before PT.
%   U_pt: the PT-ed version of U. 
%   S: cosine similarity matrix after PT.
%   S0: cosine similarity matrix before PT.
% 
% The residual is the norm of the off-diagonal elements of the cosine
% similarity matrix between the two set of vectors,  as in Eq (11) in [1].
% When comparing two pairs of vector sets, e.g., U1 and U2, the vector sets
% in U1 are more similar to each other than those in U2, if E1<E2; here the
% assumption is that both matrix sets (in U1 and U2) have the same size.
% 
% [1] https://doi.org/10.1109/OJEMB.2023.3267726
% 
% Hamid Behjat

d = inputParser;
addParameter(d, 'ApplySaturationCriteria', true);
addParameter(d, 'MaxNumberOfIterations', 20);
addParameter(d, 'SaturationParam', 1e-2); % in (0 1]
addParameter(d, 'DebugDisplay', false);
addParameter(d, 'Silent', true);
parse(d, varargin{:});
opts = d.Results;

verifyu(U);

Ns = length(U);

Ni = opts.MaxNumberOfIterations;

Nv = size(U{1},1);

Ne = size(U{1},2);

U_pt = zeros(Nv, Ne, 2);

[S0, E0] = getcossim(U);

E = zeros(Ni,1);

for iIter=1:Ni

    if iIter==1
        U_avg = U{1};
    end

    for iSet=1:Ns

        if iIter==1
            Ul = U{iSet};
        else
            Ul = U_pt(:,:,iSet);
        end

        T = getpt(U_avg, Ul);  % procrustes transform

        U_pt(:,:,iSet) = Ul * T; % apply transform

        if opts.DebugDisplay
            if iIter==1
                figure;
                axis;
                ha1 = gca;
            end
            imagesc(ha1,T);
            clim(max(clim)*[-1 1]);
            colormap(jet(101));
            colorbar;
            axis image;
            sprintf('Procrustes Transform [iteration %d/%d]',iIter,Ni);
            title(d);
            pause(0.5);
        end

    end

    if opts.DebugDisplay
        [S, E(iIter)] = getcossim(U_pt);
        if iIter==1
            figure;
            axis;
            ha2 = gca;
        end
        imagesc(ha2,S);
        clim(max(clim)*[-1 1]);
        colormap(jet(101));
        colorbar;
        axis image;
        pause(0.5);
    else
        [S, E(iIter)] = getcossim(U_pt);
    end

    if not(opts.Silent)
        prgs(iIter,Ni);
    end

    if opts.ApplySaturationCriteria
        if iIter==1
            d = E(iIter)-E0;
        else
            d = E(iIter)-E(iIter-1);
        end
        if d<0
            d = abs(d/E0);
            if d<opts.SaturationParam
                if not(opts.Silent)
                    d = 'Procrustes transform converged';
                    fprintf('\n.%s in iteration %d; .\n', d, iIter-1);
                end
                E(iIter) = 0;
                break;
            end
        end
    end

    U_avg = mean(U_pt,3);
end
end

%==========================================================================
function T = getpt(u1,u2)
d = u1'*u2;       % coordinate space bw u1 & u2
[V,~,U] = svd(d); % orthogonalize space
T = U*V';         % rotation to match u2 to u1 (u1 ~ u2*T)
end

%==========================================================================
function verifyu(U)
assert(iscell(U));
assert(all(cellfun(@(x) isequal(size(x),size(U{1})), U)));
assert(length(U)==2, 'current implementation is for only two sets.');
end

%==========================================================================
function [S, E] = getcossim(X)
if iscell(X)
    X1 = X{1};
    X2 = X{2};
else
    X1 = X(:,:,1);
    X2 = X(:,:,2);
end
C = size(X1,2);
N = repmat(vecnorm(X1)',1,C).*repmat(vecnorm(X2),C,1);
S = X1'*X2;
S = S./N; % cosine similarity
I = find(ones(C)-eye(C));
E = 0.5*sqrt(sum(S(I).^2)); % residual as in Eq (11) in [1]
end

%==========================================================================
function prgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\n Procrustes transform iteration..');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end