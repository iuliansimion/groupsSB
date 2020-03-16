#
# Read("~/Workspace/groupsSB/epi/A2/adModule.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;
L:= SimpleLieAlgebra(type,rank,Rationals);
cb0:=CanonicalBasis(L);

M:=AdjointModule(L);
extM2:=ExteriorPowerOfAlgebraModule(M,2);

ade:=function(e)
	local mb,result,v;
    
    mb:=Basis(M);

	result:=[];
	for v in mb do
		Append(result,[Coefficients(mb,e^v)]);
	od;
	#result:=TransposedMat(result);
	return result;
end;

root_group:=function(index,t)
	local e,tmp,result,i;
	
    e:=ade(cb0[index]);
    tmp:=ShallowCopy(e);
	result:=tmp^0;

	i:=1;
	while Length(Set(Concatenation(tmp)))<>1 do
		result:=result+t^i*tmp/Factorial(i);
		i:=i+1;
		tmp:=tmp*e;
	od;
	return result;
end;

ad_nil_deg:=function(e)
    local ee,tmp;
    ee:=ade(e);
    tmp:=ShallowCopy(ee);
    i:=1;
    while Length(Set(Concatenation(tmp)))<>1 do;
        tmp:=tmp*ee;
        i:=i+1;
    od;
    return i;
end;

ad_exp:=function(e,t,v)
    local result,ord_e,i,j,tmp;
    
    ord_e:=ad_nil_deg(e);

    result:=v;
    for i in [1..ord_e-1] do
        tmp:=ShallowCopy(v);
        for j in [1..i] do
            tmp:=e^tmp;
        od;
        result:=result+t^i*tmp/Factorial(i);
    od;
    return result;
end;

diag_g:=function(e,wedge_sum)
    local fam_algmod,fam_wedge,o,result,new_wedge,i,bv,exp_term;
    fam_algmod:=FamilyObj(wedge_sum);
    o:=ExtRepOfObj(wedge_sum);
    fam_wedge:=FamilyObj(o);
    o:=ExtRepOfObj(o);

    result:=[];
    for i in [1,3..Length(o)-1] do
        new_wedge:=[];
        for bv in o[i] do
            exp_term:=ad_exp(e,1,bv);
            Add(new_wedge,exp_term);
        od;
        Add(result,new_wedge);
        Add(result,o[i+1]);
    od;

    result:=ObjByExtRep(fam_wedge,result);
    result:=ConvertToNormalFormMonomialElement(result);
    result:=ObjByExtRep(fam_algmod,result);
    
    return result;
end;

exp_action_mat_g:=function(e,M)
    local mb,v,result;
    mb:=Basis(M);

    result:=[];
 	for v in mb do
		Append(result,[Coefficients(mb,e^v)]);
	od;
	#result:=TransposedMat(result);
	return result;
end;

u1:=root_group(1,1);

f:=LeftModuleHomomorphismByMatrix(Basis(M),u1,Basis(M));
test:=Image(f,Basis(M)[1]);

bbb:=Basis(extM2);
tmp:=ExtRepOfObj(ExtRepOfObj(bbb[1]));
tmp:=ShallowCopy(tmp);
tmp:=[[Image(f,tmp[1][1]),Image(f,tmp[1][2])],1];

fam:=FamilyObj(ExtRepOfObj(bbb[1]));
tmp:=ObjByExtRep(fam,tmp);
tmp:=ConvertToNormalFormMonomialElement(tmp);

fam2:=FamilyObj(bbb[1]);
tmp:=ObjByExtRep(fam2,tmp);