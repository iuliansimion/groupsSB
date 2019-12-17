
#
# Read("~/Workspace/epi/ext.gi");
#
# highest_weight:=[1,1];# A2: corresponding highest weight module is the Lie algebra
#

L:= SimpleLieAlgebra(type,rank,Rationals);
cb0:=CanonicalBasis(L);
V:= HighestWeightModule( L, highest_weight );
extA:=List([1..Dimension(V)],i->ExteriorPowerOfAlgebraModule( V, i ));
extA_basis:=List(extA,i->Basis(i));
