function U_next = SolveTriangleTomasBatch(U,G,dt,h,dim)
%Do the triangular batching job
sz = size(U);
U = shiftdim(U,dim-1);
G = shiftdim(G,dim-1);
sz = size(U);
U = reshape(U,[sz(1) sz(2)*sz(3)]);
G = reshape(G,[sz(1) sz(2)*sz(3)]);
U = U';
G = G';
for i=1:sz(2)*sz(3)
     U(i,:) = SolveTriangleTomas(squeeze(U(i,:)),squeeze(G(i,:)),dt,h);
%    U(i,:) = SolveTriangleImplicit(squeeze(U(i,:)),squeeze(G(i,:)),dt,h);
end
U = U';
G = G';
U = reshape(U,[sz(1) sz(2) sz(3)]);
G = reshape(G,[sz(1) sz(2) sz(3)]);
U = shiftdim(U,4-dim);
G = shiftdim(G,4-dim);
U_next = U;
