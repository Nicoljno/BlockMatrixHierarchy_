clear;
clc;

v = [1 2 3 4];         % Partitions
d = 4;                    % Local dimension
r = d-1;                    % SR
n = length(v);            % Number of qudits
dim = d^n;                % Dimension of the global space
id = eye(d);
for l = 0 : d-1
    comp{l+1} = id(:,l+1);
end
X = GenPauli(1,0,d);
Z = GenPauli(0,1,d);

k=[3, 2, 2];



test_results=[];

for n_test = 1:1000
    ket = RandomStateVector(dim);
    %ket = GHZState(d, n);
    %ket = sparse(dicke_state_d_dim(k));
    %ket = DickeState(length(v));
    rho_pure = ket*ket';
    vis = sdpvar(1); % Visibility

    Id = speye(dim);
    rho_mixed = vis*rho_pure + (1 - vis)*(Id/dim);
    biparts = bipartitions(v);

    Constraints = [];
    Xvars = cell(length(biparts),1);

    for i = 1:length(biparts)
        part1 = biparts{i}{1};          
        part2 = biparts{i}{2};          
        k = length(part1);             
        dim_k = d^k;
        Xvars{i} = sdpvar(dim_k, dim_k, 'hermitian', 'real');
        %Xvars{i} = diag(diag(Xvars{i}));

        Constraints = [Constraints, Xvars{i} >= 0];
        Constraints = [Constraints, eye(dim_k) - Xvars{i} >= 0];
    end

    expr_for_trace_sum = 0;
    for i = 1:length(Xvars)
        expr_for_trace_sum = expr_for_trace_sum + trace(Xvars{i});
    end
    Constraints = [Constraints, expr_for_trace_sum <= r];
    M_vars = cell(length(Xvars),1);

    dims_vec = d*ones(1,n);  % [d, d, ..., d]

    for i = 1:length(Xvars)
        part1 = biparts{i}{1};
        part2 = biparts{i}{2};
        k = length(part1);
        dim_k = d^k;
        dim_comp = d^(length(part2));
        I_comp = speye(dim_comp);
        Xi_ext = kron(Xvars{i}, I_comp);
        ordered_sites = [part1, part2];
        M_vars{i} = PermuteSystems(Xi_ext, ordered_sites, dims_vec);
    end

    M_sum = 0;
    for i = 1:length(M_vars)
        M_sum = M_sum + M_vars{i};
    end

    Constraints = [Constraints, M_sum - rho_mixed >= 0];
    Constraints = [Constraints, Id - M_sum >= 0];
    Constraints = [Constraints, 0 <= vis <= 1];

    Objective = -vis;  

    ops = sdpsettings('solver','mosek','verbose',1);
    sol = optimize(Constraints, Objective, ops);

    if sol.problem == 0
        disp('Optimization solved successfully!');
        disp(['Optimal visibility = ', num2str(value(vis))]);
    else
        disp('Something went wrong:');
        sol.info
        yalmiperror(sol.problem)
    end

    max_F=0;
    biparts = all_bipartitions(v);
    for i = 1:length(biparts)
        tmp_eig = sort(eig(PartialTrace(rho_pure, biparts{i}{1}, dims_vec)));
        tmp=0;
        for j = 1:r
            tmp=tmp+tmp_eig(end-j+1);
        end
        if(tmp>max_F)
            max_F=tmp;
        end
    end
    v_fidelity=(max_F-1/dim)/(1-1/dim);
    
    vis=value(vis);

%     yalmip clear
%     biparts = bipartitions(v);
%     vis_ =sdpvar(1);
%     rho_mixed_ = vis_*rho_pure + (1 - vis_)*(Id/dim);
%     tmp=0;
%     F=[];
%     for i=1:length(biparts)
%         part1 = biparts{i}{1};
%         part2 = biparts{i}{2};
%         ordered_sites = [part1, part2];
%         ordered_sites1 = [part2, part1];
%         sigma{i} = sdpvar(dim, dim, 'hermitian', 'complex');
%         F=[F, sigma{i}>=0];
%         tmp=tmp+sigma{i};
%         F=[F, PermuteSystems( Tensor(eye(d^length(biparts{i}{1})), PartialTrace(sigma{i}, biparts{i}{1}, dims_vec) ) , ordered_sites, dims_vec)-1/r*sigma{i}>=0 ];
%         F=[F, PermuteSystems( Tensor(eye(d^length(biparts{i}{2})), PartialTrace(sigma{i}, biparts{i}{2}, dims_vec) ) , ordered_sites1, dims_vec)-1/r*sigma{i}>=0 ];
%     end
%     F=[F, tmp==rho_mixed_];
% 
%     Objective=-vis_;
%     ops = sdpsettings('solver','mosek','verbose',1);
%     sol = optimize(F, Objective, ops);
%     v_red_map = value(vis_)
% 
%     yalmip clear
%     vis_ =sdpvar(1);
%     rho_mixed_ = vis_*rho_pure + (1 - vis_)*(Id/dim);
%     F=[PartialTranspose(rho_mixed_,1)>=0];
%     Objective=-vis_;
%     ops = sdpsettings('solver','mosek','verbose',1);
%     sol = optimize(F, Objective, ops);
%     v_ppt = value(vis_)

    % yalmip clear
    % Z = sdpvar(d,d,'hermitian','complex');
    % vis_ =sdpvar(1);
    % rho_mixed_ = vis_*rho_pure + (1 - vis_)*(Id/dim);
    % F=[Tensor(Z, eye(d))-rho_mixed_>=0];
    % F=[F, trace(Z)/d<=r/d];
    % Objective=-vis_;
    % ops = sdpsettings('solver','mosek','verbose',1);
    % sol = optimize(F, Objective, ops);
    % v_fidelity_2 = value(vis_);

    
%     yalmip clear 
%     vis_ =sdpvar(1);
%     rho_a=sdpvar(d,d,'hermitian','complex')
%     Pi1=sdpvar(d,d,'hermitian','complex')
%     rho_mixed = vis_*rho_pure + (1 - vis_)*(Id/dim);
%     
%     F=[rho_a>=0, rho_b>=0];
%     F=[F; Pi1>=0, trace(Pi1)==1];
%     F=[F, rho_a==PartialTrace(rho_mixed, 2)];
%     F=[F, trace(rho_a)==1];

    
%     Gamma=[[Tensor(Pi1,eye(d)) Tensor(rho_a, eye(d)) rho_mixed]; ... 
%         [Tensor(rho_a, eye(d)) Tensor(rho_a, eye(d)) rho_mixed]; ...
%         [rho_mixed rho_mixed rho_mixed];];
% 
%     
%     F=[F, Gamma>=0];
%     Objective=-vis_;
%     ops = sdpsettings('solver','mosek','verbose',1);
%     sol = optimize(F, Objective, ops);
%     v_fidelity_3 = value(vis_);

    test_results=[test_results;  vis v_fidelity]

end


%%  Function bipartitions(v): returns the minimal non-repeating bipartitions of v
function parts = bipartitions(v)
    n = length(v);
    parts = {};
    % Only go up to 2^(n-1)-1 to skip empty set and full set.
    for mask = 1:(2^(n-1)-1)
        part1 = [];
        part2 = [];
        for i = 1:n
            if bitget(mask, i) % check if i-th bit is set
                part1 = [part1 v(i)];
            else
                part2 = [part2 v(i)];
            end
        end
        % Sort the two subsets by length (just to be consistent with the Julia code)
        if length(part1) <= length(part2)
            parts{end+1} = {part1, part2};
        else
            parts{end+1} = {part2, part1};
        end
    end
end

%%  Function all_bipartitions(v): returns all the bipartitions of v
function parts = all_bipartitions(v)
    n = length(v);
    parts = {};
    % Only go up to 2^(n-1)-1 to skip empty set and full set.
    for mask = 1:(2^(n-1)-1)
        part1 = [];
        part2 = [];
        for i = 1:n
            if bitget(mask, i) % check if i-th bit is set
                part1 = [part1 v(i)];
            else
                part2 = [part2 v(i)];
            end
        end
        parts{end+1} = {part1, part2};
        parts{end+1} = {part2, part1};
    end
end

%% function DickeState(k)
function DickeState = dicke_state_d_dim(k)
% DICKE_STATE_D_DIM  Construct a d-dimensional Dicke state for n qudits.
%
%   DickeState = dicke_state_d_dim(k)
%
%   k is a 1xd vector [k0, k1, ..., k_{d-1}] whose entries sum to n.
%   The result DickeState is a column vector in C^(d^n) representing
%   the Dicke state
%
%       |D> = (1 / sqrt(#S)) * sum_{w in S} |w>,
%
%   where S is the set of all distinct permutations of the multiset
%   determined by k.
%
%   Example:
%       % Suppose d=3 and k=[2,1,1] (so n=4).
%       % That corresponds to the multiset {0,0,1,2}.
%       D = dicke_state_d_dim([2, 1, 1])

    % Make sure k sums to n.
    n = sum(k);
    d = length(k);
    
    % Create the base multiset array. For example, if k = [2,1,1] and d=3,
    % arr = [0, 0, 1, 2].
    arr = [];
    for val = 1:d
        arr = [arr, (val-1)*ones(1, k(val))];
    end

    % Find all distinct permutations of this multiset.
    % Each row of P is one distinct permutation.
    P = unique(perms(arr), 'rows');

    % The total dimension is d^n. We build the state as a d^n x 1 vector.
    dim = d^n;
    DickeState = zeros(dim,1);

    % For each distinct permutation, compute its linear index in the
    % standard computational basis and increment the corresponding entry.
    % (Recall that MATLAB uses 1-based indexing.)
    for iPerm = 1:size(P,1)
        idx = 1;  % Start with 1 for 1-based indexing
        for position = 1:n
            idx = (idx - 1)*d + P(iPerm, position) + 1;
        end
        DickeState(idx) = DickeState(idx) + 1;
    end

    % Normalize by the square root of the number of distinct permutations:
    % #S = n! / (k_0! k_1! ... k_{d-1}!).
    % We can also just use size(P,1) for #S.
    DickeState = DickeState / sqrt(size(P,1));
end
