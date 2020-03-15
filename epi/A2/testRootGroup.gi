
#
# Read("~/Workspace/groupsSB/epi/A2/testRootGroup.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

L:= SimpleLieAlgebra(type,rank,Rationals);
cb0:=CanonicalBasis(L);
V:= HighestWeightModule( L, highest_weight );
extA:=List([1..Dimension(V)],i->ExteriorPowerOfAlgebraModule( V, i ));
extA_basis:=List(extA,i->Basis(i));

#
# change basis order of first exteriror power to match chevalley basis order
#
b:=extA_basis[1];
#SetBasis(extA[1],bbb);
#extA_basis[1]:=Basis(extA[1]);
#extA_basis[1]:=
#b:=[
#    b[7],
#    b[6],
#    b[8],   # alpha + beta
#    b[3],
#    b[2],
#    b[1],
#    b[5],
#    b[4],
#    ];

#plist:=[7,6,8,3,2,1,5,4];
plist:=[3,2,1,6,7,8,4,5];
bbb:=b{plist};
apply_plist:=function(mat)
    local result;
    result:=List(mat,i->i{plist});
    result:=TransposedMat(result);
    result:=List(result,i->i{plist});
    result:=TransposedMat(result);
    return result;
end;

semne:=DiagonalMat([1,-1,1,-1, 1,1,-1,1]);
apply_signs:=function(mat)
    local result;
    result := semne*mat*semne;
    return result;
end;

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

ade:=function(e)
	local result,v;
	result:=[];
	for v in cb0 do
		Append(result,[Coefficients(cb0,e*v)]);
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

ad_root_group:=function(index,t)
	local ee,tmp,result,i;
	ee:=ade(cb0[index]);
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


#
# Test:
#

# ade1
yyy:=ext1e(cb0[1]);
yyy:=apply_plist(yyy);
yyy:=apply_signs(yyy);
xxx:=ade0(cb0[1]);
Display(yyy);
Display(xxx);
Display(xxx=yyy);

# root group 1
yyy:=ext1_root_group(1,1);
yyy:=apply_plist(yyy);
yyy:=apply_signs(yyy);
xxx:=ad_root_group(1,1);
Display(yyy);
Display(xxx);
Display(xxx=yyy);

# ade2 --------------------> Distinct
yyy:=ext1e(cb0[2]);
yyy:=apply_plist(yyy);
yyy:=apply_signs(yyy);
xxx:=ade0(cb0[2]);
Display(yyy);
Display(xxx);
Display(xxx=yyy);

# root group 2 --------------------> Distinct
yyy:=ext1_root_group(2,1);
yyy:=apply_plist(yyy);
yyy:=apply_signs(yyy);
xxx:=ad_root_group(2,1);
Display(yyy);
Display(xxx);
Display(xxx=yyy);
