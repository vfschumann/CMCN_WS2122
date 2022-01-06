function out_mat = coupling(model)
    % function computing the coupling type between all reactions
    % in the model
    % Input: reduced model, all irreversible
    % Output: matrix with codes about coupling state:
    % 0 - not coup, 1 - fully coup, 2- partially coup, 3-direct. coup i->j,
    % 4 - direct. coup j -> i
    
    % TODO: add check for input format assumptions
    
    suff_large_num = 10^8;
    
    out_mat = ones(size(model.rxns, 1), size(model.rxns, 1)) * -1;
    
    Aeq = model.S;
    Aeq = [Aeq, zeros(size(model.mets))];

    beq = zeros(size(Aeq, 1), 1);

    A = eye(size(model.rxns, 1));        
    A = [A, -model.ub];

    b = zeros(size(A, 1), 1);
        
    for rxn_idx_from = 1:size(model.rxns, 1)
        for rxn_idx_to = 1:size(model.rxns, 1)
            
            % before doing anything check for associative properties
            % and whether they allow us to determine the coupling of this 
            % pair if the reverse was already calculated
            assoc_coup = [1 2];
            reverse_coup = out_mat(rxn_idx_to, rxn_idx_from);
            if reverse_coup ~= -1
                if ismember(reverse_coup, assoc_coup)
                    out_mat(rxn_idx_from, rxn_idx_to) = reverse_coup;
                end
            end

            lb = zeros(size(model.rxns, 1)+ 1, 1);

            ub = ones(size(model.rxns, 1)+ 1, 1) * suff_large_num;

            f = model.c;
            f = [f; 0];

            lb(rxn_idx_to) = 1;
            ub(rxn_idx_to) = 1;

            f(rxn_idx_from) = 1;
   
            lin_res_min = linprog( f, A, b, Aeq, beq, lb, ub);
            lin_res_max = linprog(-f, A, b, Aeq, beq, lb, ub);  
            
            min = lin_res_min(rxn_idx_from);
            max = lin_res_max(rxn_idx_from);
            
            % determine type from min and max
            if min == 0
                if max == suff_large_num
                    % not coupled
                    coup = 0;
                elseif max > 0
                    % directional i -> j
                    coup = 3;
                end
            elseif min > 0
                if max == suff_large_num
                    % directional j -> i
                    coup = 4;
                elseif max == min
                    % full
                    coup = 1;                
                elseif max > 0
                    % partial
                    coup = 2;
                end
            end
            
            out_mat(rxn_idx_from, rxn_idx_to) = coup;
            
        end
    end   
end