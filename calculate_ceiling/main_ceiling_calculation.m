% this is from your spreadsheet: cols are the 14 clades; rows are the weeks

%UK
%sample_counts = [23075 3033 3601 1995 8 41 3956 51395 804 1386 223 153 20 11];
%USA
sample_counts = [7132 227 87 402 25 39 1407 12192 1151 32860 18 74 1040 445];



% this is from your Figure 1B. Order of cols correspond with those in
% 'sample_counts' above

clade_defining_SNP_diffs = [0 10 2 7 12 10 5 3 4 5 14 9 12 16;
    10 0 12 3 6 4 9 7 8 9 8 5 8 12;
    2 12 0 9 14 12 7 5 6 7 16 11 14 18;
    7 3 9 0 5 3 6 4 5 6 7 2 5 9;
    12 6 14 5 0 6 11 9 10 11 12 7 10 14;
    10 4 12 3 6 0 9 7 8 9 10 5 8 12;
    5 9 7 6 11 9 0 2 1 2 13 8 11 15;
    3 7 5 4 9 7 2 0 1 2 11 6 9 13;
    4 8 6 5 10 8 1 1 0 1 12 7 10 14;
    5 9 7 6 11 9 2 2 1 0 13 8 11 15;
    14 8 16 7 12 10 13 11 12 13 0 5 8 12;
    9 5 11 2 7 5 8 6 7 8 5 0 3 7;
    12 8 14 5 10 8 11 9 10 11 8 3 0 10;
    16 12 18 9 14 12 15 13 14 15 12 7 10 0];

% the clade-defining SNP threshold to call a recombinant:
% if SNP_threshold = 3, for example, a sequence needs to have at least
% three clade-defining SNPs from one parent and at least three clade-defining SNPs from the other parent at the sites where they differ 
SNP_threshold = 3; 
n_SNPs = 37;    % not needed in code anywhere

[n_weeks, n_clades] = size(sample_counts);

%recomb_sample_counts = [46]; %UK
recomb_sample_counts = [48]; %USA

sample_size_by_week = sum(sample_counts, 2); % just calculates how many samples per week

clade_freqs = sample_counts./sample_size_by_week; % calculates the frequencies of each clade by week

prob_observing_recomb = GetProbabilityObservingRecombinants(n_clades, clade_defining_SNP_diffs, SNP_threshold); % this function returns a 14 x 14 matrix with the probability of a recombinant being observed if it happened between clade i and clade j
writematrix(prob_observing_recomb,"prob_obs_cdSNP-3.txt")

% for each week, find highest proportion of the viral population that could be recombinant, where the 95% confidence interval lower bound still contains the observed number of recombinants
for week_num = 1:n_weeks
    this_week_freq = clade_freqs(week_num, :);
    highest_recomb_prop(week_num) = GetHighestRecombProportion(this_week_freq, prob_observing_recomb, n_clades, sample_size_by_week(week_num), recomb_sample_counts(week_num));
end
writematrix(highest_recomb_prop,"recomb_prop_cdSNP-3.USA.txt")


% redo-ing the analysis with a SNP threshold of 2 instead of 3:

% try a threshold of 2
SNP_threshold = 2; % to call a recombinant, it has to have at least two SNPs uniquely present from the minor parent
n_SNPs = 37;

%recomb_sample_counts = [100]; %UK
recomb_sample_counts = [595]; %USA

sample_size_by_week = sum(sample_counts, 2);

clade_freqs = sample_counts./sample_size_by_week; 

prob_observing_recomb = GetProbabilityObservingRecombinants(n_clades, clade_defining_SNP_diffs, SNP_threshold);
writematrix(prob_observing_recomb,"prob_obs_cdSNP-2.txt")

for week_num = 1:n_weeks
    
    %find highest proportion_recombinant whose 95% confidence interval still contains the observed number of recombinants
    this_week_freq = clade_freqs(week_num, :);
    highest_recomb_prop(week_num) = GetHighestRecombProportion(this_week_freq, prob_observing_recomb, n_clades, sample_size_by_week(week_num), recomb_sample_counts(week_num));
    [week_num highest_recomb_prop(week_num)]
end
writematrix(highest_recomb_prop,"recomb_prop_cdSNP-2.USA.txt")

