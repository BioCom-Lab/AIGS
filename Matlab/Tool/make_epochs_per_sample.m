function [result] = make_epochs_per_sample(weight)
result = -1*ones(size(weight));
result(weight>0) = max(weight,[],'all')./weight(weight>0);
% n_samples = n_epochs*(weight/max(weight,[],'all'));
% result(n_samples>0) = n_epochs./n_samples(n_samples>0);
end