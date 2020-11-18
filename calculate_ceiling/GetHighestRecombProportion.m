function highest_recomb_prop = GetHighestRecombProportion(this_week_freq, prob_observing_recomb, n_clades, tot_samples, n_recomb_observed)

%this_week_freq
%pause

if (tot_samples == 0) || isnan(this_week_freq(1))
    highest_recomb_prop = NaN; return;
end

% start off at ceiling on recombinant proportion at 0 and move up:
prop_recombinant = 0;
prop_incr = 0.00001;

while 1
    overall_prob = 0;
    for i = 1:n_clades
        for j = 1:n_clades
            % calculate probability that parental strains are i and j:
            prob_parental_strains = this_week_freq(i)*this_week_freq(j);
            % overall prob keeps track of the overall probability of
            % observing a recombinant in the focal week, given that the
            % proportion of the viral population that is recombinant is prop_recombinant
            overall_prob = overall_prob + prop_recombinant*prob_parental_strains*prob_observing_recomb(i,j);
        end
    end
    % X contains the 95% confidence interval on the number of recombinants that are observed under the specified number of samples and given the overall probability of observing a recombinant
    X = binoinv([0.025 0.975],tot_samples,overall_prob);
    %[prop_recombinant X(1) n_recomb_observed tot_samples]
    % if the lower bound excedes the number of recombinants observed, prop_recombinant is too high - go with last one in the loop 
    if X(1) > n_recomb_observed
        break;
    end
    if prop_recombinant > 1
        highest_recomb_prop = 1;
        break;
    end
    highest_recomb_prop = prop_recombinant;
    prop_recombinant = prop_recombinant + prop_incr;
end

