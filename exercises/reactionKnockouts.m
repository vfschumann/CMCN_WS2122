function vec = reactionKnockouts(cycle, solved)
    
    limit = length(cycle.rxns);
    arr = zeros(limit,1);
    vec = arr;

    z_opt = solved.f;
    
    for q=1:limit

        id = cycle.rxns(q);
           
        upd_solution = changeRxnBounds(cycle, id, [0],'b');
        opt_solution = optimizeCbModel(upd_solution);     
        vec(q) = opt_solution.f / z_opt;
    end     
    return
end