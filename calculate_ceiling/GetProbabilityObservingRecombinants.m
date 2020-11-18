function prob_observing_recomb = GetProbabilityObservingRecombinants(n_clades, clade_defining_SNP_diffs, SNP_threshold)

prob_observing_recomb = zeros(n_clades);

for i = 1:n_clades
    for j = i:n_clades
        SNP_diff = clade_defining_SNP_diffs(i,j);
        if SNP_diff < 2*SNP_threshold
            prob_observing_recomb(i,j) = 0;
        else
            % figure out the probability mass of recombinants that have at
            % least SNP_threshold unique SNPs from the minor parent; sum up
            % this mass to give the probability of calling the recombinant sequence a recombinant if the parental viruses were as specified
            y = binopdf([SNP_threshold:(SNP_diff-SNP_threshold)],SNP_diff,0.5);
            prob_observing_recomb(i,j) = sum(y); y = [];
        end
        prob_observing_recomb(j,i) = prob_observing_recomb(i,j);
    end
end
