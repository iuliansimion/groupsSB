#
# Read("~/Workspace/groupsSB/epi/A2/actionFromExt1.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;
highest_weight:=[1,1];# corresponding highest weight module is the Lie algebra

L:= SimpleLieAlgebra(type,rank,Rationals);
cb0:=CanonicalBasis(L);
V:= HighestWeightModule( L, highest_weight );
extA:=List([1..Dimension(V)],i->ExteriorPowerOfAlgebraModule( V, i ));
extA_basis:=List(extA,i->Basis(i));

extA1:=extA[1];
extA1b:=extA_basis[1];

#extA:=List([1..Dimension(V)],i->ExteriorPower( i ,V ));

ext1e:=function(e)
	local extb,result,v;
    extb:=extA_basis[1];

	result:=[];
	for v in extb do
		Append(result,[Coefficients(extb,e^v)]);
	od;
	result:=TransposedMat(result);
	return result;
end;

ext1_root_group:=function(index,t)
	local ee,tmp,result,i;
	ee:=ext1e(cb0[index]);
	tmp:=ee;
	result:=tmp^0;
	i:=1;
	while Length(Set(Concatenation(tmp)))<>1 do
		result:=result+t^i*tmp/Factorial(i);
		i:=i+1;
		tmp:=tmp*ee;
	od;
	return result;
end;

ext1_coeff_comb:=function(coeffs)
    local result,extb;
    extb:=extA_basis[1];
    result:=List([1..Length(extb)],i->coeffs[i]*extb[i]);
    return Sum(result);
end;

ext1_root_group_basis_img:=function(index,t)
    local rg,result;
    rg:=ext1_root_group(index,t);
    rg:=TransposedMat(rg);
    result:=List(rg,l->ext1_coeff_comb(l));
    return result;
end;


# ObjByExtRep(FamilyObj(bbb[1]),ExtRepOfObj(bbb[1]));

u1Bimg:=ext1_root_group_basis_img(1,1);

qqq:=ExtRepOfObj(extA_basis[2][1]);
IsWedgeElement(qqq);
ConvertToNormalFormMonomialElement(qqq);

n:=2;
fam:=FamilyObj(qqq);
combs:= Combinations( [1..Dimension(L)], n );
gens:= List( combs, x -> Basis(L){x} );
F:=LeftActingDomain(L);
gens2:= List( gens, x -> ObjByExtRep( fam, [ x , One(F) ] ) );
    
for i in [1..Length(gens)] do
    gens2[i]![2]:= true;
od;

test:=[2,2];
test:=Basis(L){test};
test:=ObjByExtRep( fam, [ test , One(F) ] );
test![2]:=true;
test:=ConvertToNormalFormMonomialElement(test);


M:=AdjointModule(L);
extA2:=ExteriorPowerOfAlgebraModule(M,2);


#
#
#
extA1:=extA[1];
bv:=BasisVectors(Basis(extA1));
bb:=ext1_root_group_basis_img(1,1);
B:= Algebra( Rationals, bb, "basis" );
#AlgebraHomomorphismByImages( extA1, extA1, bv, bb )

LeftModuleGeneralMappingByImages( extA1, extA1, bv, bb );
LeftModuleHomomorphismByImages( extA1, extA1, bv, bb );