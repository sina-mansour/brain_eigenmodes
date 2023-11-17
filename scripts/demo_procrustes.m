%
% 1. The script requires two functions:
% .https://github.com/aitchbi/matlab-utils/blob/main/hb_procrustes.m
% .https://github.com/NSBLab/BrainEigenmodes/blob/main/functions_matlab/bluewhitered.m
% Copy & addpath. 
%
% 2. Eigenmodes are loaded as NxM matrices, where M is the number or
% eigenmodes; U1.mat, etc. Four sample sets of eigenmodes can be downloaded
% from here to test the script:
% https://www.dropbox.com/scl/fo/5brk89ypsdywuhf98f5k2/h?rlkey=aevjvak5lss6odfvyllsk6ync&dl=0
%
% 3. See NOTE 1 & NOTE 2 for the two settings. 

clc
clear
close all


load('eigenmodes.mat')

% f_U{1} = 'U1.mat';
% f_U{2} = 'U2.mat';
% f_U{3} = 'U3.mat';
% f_U{4} = 'U4.mat';

f_U{1} = geometry_eigenmodes;
f_U{2} = edr_eigenmodes;
f_U{3} = pang_connectome_eigenmodes;
f_U{4} = our_connectome_eigenmodes;


%...extend list as needed. 

WhichPairs = [1 2 3]; % NOTE 1

% HowManyModes = 10:10:100; %10:10:100; % NOTE 2

HowManyModes = 2:1:200; %10:10:100; % NOTE 2
% HowManyModes = 100; %10:10:100; % NOTE 2

Pairs = {
    '01-vs-02'
    '01-vs-03'
    '01-vs-04'
    %...extend pairs as needed.
    };

Pairs = Pairs(WhichPairs);

% NOTE 1
% Which pair of eigenmodes to compare, e.g., 1 mean the first set in Pairs.
% And within the Pairs list, 'a-vs-b' means comparing f_U{a} and f_U{b}; a
% and b \in 00-99.
%
% NOTE 2
% How many of initial eigenmodes to compare, e.g., [A B C] means to compare
% three sets of eigenmodes: the first A, the first B, and the first C. This
% is done in each of the pair sets. If you set this variable to a single
% value (e.g., 100), the script generates the cosine similarty matrix
% plots, otherwise the "residual vs number of eigenmodes" plot is
% generated.

N_mode_sets  = length(HowManyModes);

N_pairs      = length(WhichPairs);

for iPair=1:N_pairs

    ThePair = Pairs{WhichPairs(iPair)};

    E0_sets = zeros(N_mode_sets,1);
    E_sets  = zeros(N_mode_sets,1);
    P_sets  = zeros(N_mode_sets,1);

    % [U1, U2] = loadmodes(ThePair, f_U);

    [U1, U2] = myloadmodes(ThePair, f_U);
    % U1 = f_U(str2num(ThePair(1:2)));
    % U1 = U1{1};
    % U2 = f_U(str2num(ThePair(7:8)));
    % U2 = U2{1};

    for iSet=1:N_mode_sets

        Ne = HowManyModes(iSet);

        U = cell(1,2);
        U{1} = U1(:,1:Ne);
        U{2} = U2(:,1:Ne);

        if N_mode_sets==1

            [E, E0, ~, S, S0] = hb_procrustes(U);

            d = find(E);
            E = E(d(end));

            disp(' ');
            fprintf('\n.Residual, before PT: %0.04f', E0);
            fprintf('\n.Residual,  after PT: %0.04f', E);
            fprintf('\n.Percentage of residual left: %0.04f%%\n', E/E0*100);
            disp(' ');

            hf = figure;
            set(hf, 'Position', [10 10 1000 400])
            subplot(121);
            imagesc(S0);
            axis image;
            title(['cosine similarity before PT (', ThePair, ')' ]);
            xlabel('eigenmodes'); ylabel('eigenmodes');
            set(gca, 'FontSize',12);
            clim(max(clim)*[-1 1]); colormap(bluewhitered(101)); colorbar;
            subplot(122);
            imagesc(S);
            axis image;
            title(['cosine similarity after PT (', ThePair, ')' ]);
            xlabel('eigenmodes'); ylabel('eigenmodes');
            set(gca, 'FontSize',12);
            clim(max(clim)*[-1 1]); colormap(bluewhitered(101)); colorbar;

        else

            [E, E0, ~, S, S0] = hb_procrustes(U);

            d = norm(S, 'fro') - norm(S0, 'fro');

            assert(abs(d) < 1e-6)

            d = find(E);
            E = E(d(end));

            E0_sets(iSet) = E0 * 2 / norm(S0, 'fro') * 100;
            E_sets(iSet)  = E * 2 / norm(S, 'fro') * 100;;
            P_sets(iSet) = E/E0*100;

            if iSet==N_mode_sets
                if 0 %1
                    % individual plots for each pair
                    figure;
                    subplot(311);
                    plot(HowManyModes,E0_sets, '-ok', 'MarkerFaceColor','k');
                    subplot1(ThePair);
                    subplot(312);
                    plot(HowManyModes,E_sets, '-ok', 'MarkerFaceColor','k');
                    subplot2(ThePair);
                    subplot(313);
                    plot(HowManyModes,P_sets, '-ok', 'MarkerFaceColor','k');
                    subplot3(ThePair);
                end
            end
            fprintf('\n.Pair %d/%d, set %02d/%02d done.',...
                iPair, N_pairs, iSet, N_mode_sets);
        end
    end

    if N_pairs>1 && N_mode_sets>1

        if iPair==1
            E0_pairs = cell(N_pairs,1);
            E_pairs  = cell(N_pairs,1);
            P_pairs  = cell(N_pairs,1);
            settings = cell(N_pairs,1);
        end
        E0_pairs{iPair} = E0_sets;
        E_pairs{iPair}  = E_sets;
        P_pairs{iPair}  = P_sets;
        settings{iPair} = ThePair;

        if iPair==N_pairs
            figure;
            subplot(311);
            hold on;
            for k=1:N_pairs
                plot(HowManyModes,E0_pairs{k}, '-o');
            end
            legend(settings,'Location','southeast');
            subplot1;
            subplot(312);
            hold on;
            for k=1:N_pairs
                plot(HowManyModes,E_pairs{k}, '-o');
            end
            legend(settings,'Location','southeast');
            subplot2;
            subplot(313);
            hold on;
            for k=1:N_pairs
                plot(HowManyModes,P_pairs{k}, '-o');
            end
            legend(settings,'Location','southeast');
            subplot3;
        end
    end
end

% save outputs to file
save('procrustes_results.mat', "HowManyModes", "E0_pairs", "E_pairs", "P_pairs", "Pairs")

fprintf('\n.Script done.')



function [U1, U2] = myloadmodes(t, f_U)
    assert(length(t)==8);
    assert(strcmp(t(3:6),'-vs-'));
    U1 = f_U(str2num(t(1:2)));
    U1 = U1{1};
    U2 = f_U(str2num(t(7:8)));
    U2 = U2{1};
end


function [U1, U2] = loadmodes(t,f_U)
assert(length(t)==8);
assert(strcmp(t(3:6),'-vs-'));
a = str2double(t(1:2));
b = str2double(t(7:8));
d = load(f_U{a});
U1 = d.U;
d = load(f_U{b});
U2 = d.U;
end

function subplot1(t)
grid on;
xlabel('number of modes');
ylabel('%');
d = '[% of Frobenius norm]';
if exist('t', 'var')
    title(['pre procrustes tranform residual ', d, ' (', t, ')']);
else
    title('pre procrustes tranform residual ', d);
end
end

function subplot2(t)
grid on;
xlabel('number of modes');
ylabel('%');
d = '[% of Frobenius norm]';
if exist('t', 'var')
    title(['post procrustes tranform residual ', d, ' (', t, ')']);
else
    title('post procrustes tranform residual ', d);
end
end

function subplot3(t)
grid on;
xlabel('number of modes');
ylabel('%');
if exist('t', 'var')
    title(['% of residual remaining after procrustes transform (', t, ')']);
else
    title('% of residual remaining after procrustes transform');
end
end


