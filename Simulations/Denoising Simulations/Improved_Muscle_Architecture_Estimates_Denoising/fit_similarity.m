function [S, Rcs, D, Dpos] = fit_similarity(fiber_mm_a, fiber_mm_b, C)
%FIT_SIMILARITY Compute similarity between two fibers
%   fiber_mm_a: N-by-3 matrix of fiber A coordinates in *MILLIMETER*
%   fiber_mm_b: N-by-3 matrix of fiber B coordinates in *MILLIMETER*
%   C: Weighting factor in similarity equation
%
% OUTPUT ARGUMENTS
%   S: Similarity coefficient
%   Rcs: Corresponding segment ratio
%   D: Mean Euclidean Distance
%   Dpos: Positional difference within correponding segment (a - b)
%
% *NOTE*: Fiber A and fiber B must be sampled at identical step size and
% share identical starting point, otherwise the result may be erroneous.
%
% Reference:
%   Ding et al. (2001) Case Study: Reconstruction, Visualization And
%   Quantification Of Neuronal Fiber Pathways

% Vector length
vlen_a = size(fiber_mm_a, 1);
vlen_b = size(fiber_mm_b, 1);

% Curve length
clen_a = get_full_len(fiber_mm_a);
clen_b = get_full_len(fiber_mm_b);

if vlen_a < vlen_b
    Rcs = clen_a / clen_b;
    fiber_mm_b = fiber_mm_b(1 : vlen_a, : );
else
    Rcs = clen_b / clen_a;
    fiber_mm_a = fiber_mm_a(1 : vlen_b, : );
end

% Mean Euclidean Distance
Dpos = fiber_mm_a - fiber_mm_b;
D = mean(sqrt(sum(Dpos .^ 2, 2)));

S = Rcs .* exp(-D / C);
end