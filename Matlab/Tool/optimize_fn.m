function [head_embedding,epochs_per_sample,epochs_per_negative_sample,epoch_of_next_sample,epoch_of_next_negative_sample] = optimize_fn(head_embedding,tail_embedding,head,tail,n_vertices,epochs_per_sample,a,b,alpha,epochs_per_negative_sample,epoch_of_next_negative_sample,epoch_of_next_sample,n,A_att,A_rep)
%仿照layouts中_optimize_layout_euclidean_single_epoch()

for i = 1:size(epochs_per_sample,1)
    if epoch_of_next_sample(i)<=n
        j = head(i);
        k = tail(i);
        current = head_embedding(j,:);
        other = tail_embedding(k,:);
        dist_squared = norm(current-other);
        dist_squared = dist_squared^2;
        if dist_squared>0
            grad_coeff = -2*a*b*dist_squared^(b-1);
            grad_coeff = grad_coeff/(a*dist_squared^b+1);
        else
            grad_coeff = 0;
        end
        for d = 1:2
            grad_d = clip(grad_coeff*(current(d)-other(d)));
            current(d) = current(d)+alpha*grad_d;
        end
        epoch_of_next_sample(i) = epoch_of_next_sample(i)+epochs_per_sample(i);
        epochs_per_sample(i) = epochs_per_sample(i)+epochs_per_sample(i);
        n_neg_samples = floor((n-epoch_of_next_negative_sample(i))/epochs_per_negative_sample(i));
        for p = 1:n_neg_samples
            k = randi(n_vertices);
            if k == j
                break
            end
            other = tail_embedding(k,:);
            dist_squared = norm(current-other);
            dist_squared = dist_squared^2;
            if dist_squared>0
                grad_coeff = 2*b;
                grad_coeff = grad_coeff/(0.001+dist_squared)/(a*dist_squared^b+1);
            elseif j==k
                continue
            else
                grad_coeff = 0;
            end
            for d = 1:2
                if grad_coeff>0
                    grad_d = clip(grad_coeff*(current(d)-other(d)));
                else
                    grad_d = 4;
                end
                current(d) = current(d)+alpha*grad_d;
            end
            head_embedding(j,:) = current;
            epoch_of_next_negative_sample(i) = epoch_of_next_negative_sample(i)+n_neg_samples*epochs_per_negative_sample(i);
            epochs_per_negative_sample(i) = epochs_per_negative_sample(i)+n_neg_samples*epochs_per_negative_sample(i);
        end
    end
end