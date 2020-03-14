
#
# Read("~/Workspace/groupsSB/epi/ext.gi");
#
# highest_weight:=[1,1];# A2: corresponding highest weight module is the Lie algebra
#

L:= SimpleLieAlgebra(type,rank,Rationals);
cb0:=CanonicalBasis(L);
V:= HighestWeightModule( L, highest_weight );
extA:=List([1..Dimension(V)],i->ExteriorPowerOfAlgebraModule( V, i ));
extA_basis:=List(extA,i->Basis(i));

ext_e:=function(k,e)
	local W,extb,result,v;
    W:=extA[k];
    extb:=extA_basis[k];

	result:=[];
	for v in extb do
		Append(result,[Coefficients(extb,e^v)]);
	od;
	result:=TransposedMat(result);
	return result;
end;

ade0:=function(e)
	local result,v;
	result:=[];
	for v in cb0 do
		Append(result,[Coefficients(cb0,e*v)]);
	od;
	result:=TransposedMat(result);
	return result;
end;


ext_root_group:=function(k,index,t)
	local ee,tmp,result,i;
	ee:=ext_e(k,cb0[index]);
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

ext_root_nil_degree:=function(k,index)
	local ee,tmp,i;
	ee:=ext_e(k,cb0[index]);
	tmp:=ee;
	i:=1;
	while Length(Set(Concatenation(tmp)))<>1 do
		i:=i+1;
		tmp:=tmp*ee;
	od;
	return i;
end;

ext_lin_comb:=function(k,coeffs)
	local extAb,result;
	extAb:=extA_basis[k];
	result:=List([1..Length(extAb)],i->coeffs[i]*extAb[i]);
	return Sum(result);
end;