
#
# Read("~/Workspace/epi/group.gi");
#
#
# needs:
# type:="A";
# rank:=2;
# nr_pos_roots:=3;
#


ZZ:=Integers;
avarnames:=List([1..100],i->Concatenation("a_{",String(i),"}"));
bvarnames:=List([1..100],i->Concatenation("b_{",String(i),"}"));
cvarnames:=List([1..100],i->Concatenation("c_{",String(i),"}"));
xvarnames:=List([1..100],i->Concatenation("x_{",String(i),"}"));
varnames:=Concatenation(avarnames,bvarnames,cvarnames,xvarnames);
APR:=PolynomialRing(ZZ,varnames);
vars:=IndeterminatesOfPolynomialRing(APR);
xvars:=vars{[301..400]};

sla:=SimpleLieAlgebraTypeA_G(type,rank,APR);

cb:=CanonicalBasis(sla);

e:=cb[1];

id_mat:=DiagonalMat(List([1..2*nr_pos_roots+rank],i->1));


ade:=function(e)
	local result,v;
	result:=[];
	for v in cb do
		Append(result,[Coefficients(cb,e*v)]);
	od;
	result:=TransposedMat(result);
	return result;
end;


root_group:=function(index,t)
	local ee,tmp,result,i;
	ee:=ade(cb[index]);
	tmp:=ee;
	result:=One(APR)*tmp^0;
	i:=1;
	while Length(Set(Concatenation(tmp)))<>1 do
		result:=result+t^i*tmp/Factorial(i);
		i:=i+1;
		tmp:=tmp*ee;
	od;
	return result;
end;
#u1a1:=root_group(1,vars[1]);


pos_root_groups:=function(start_a_index)
	return List([1..nr_pos_roots],i->root_group(i,vars[start_a_index+i]));
end;

generic_U:=function(start_a_index)
	local Uas;
	Uas:=pos_root_groups(start_a_index);
	return Product(Uas);
end;

Ua:=generic_U(0);
Ub:=generic_U(10);
Uc:=generic_U(20);

#
#
#

evaluate_U:=function(u,vals)
	local i,j,result,v;
	result := [];
	for i in [1..Length(u)] do
		Append(result,[[1..Length(u)]]);
		for j in [1..Length(u)] do
			result[i][j]:=u[i][j];
		od;
	od;
	Print(result);
	for i in [1..Length(u)] do
		for j in [1..Length(u)] do
			for v in vals do
				result[i][j]:=One(APR)*Value(One(APR)*result[i][j],v[1],v[2]);
			od;
		od;
	od;
	return result;
end;



evaluate_U:=function(u,vals)
	local i,j,result,v;
	result := [];
	for i in [1..Length(u)] do
		Append(result,[[1..Length(u)]]);
		for j in [1..Length(u)] do
			result[i][j]:=u[i][j];
		od;
	od;
	Print(result);
	for i in [1..Length(u)] do
		for j in [1..Length(u)] do
			for v in vals do
				result[i][j]:=One(APR)*Value(One(APR)*result[i][j],v[1],v[2]);
			od;
		od;
	od;
	return result;
end;

evaluate_rels:=function(rels,vals)
	local i,result,v;
	result :=List([1..Length(rels)],i->rels[i]);
	for i in [1..Length(rels)] do
		for v in vals do
			#Print(Length(rels),": ",rels[i],"\n");
			result[i]:=One(APR)*Value(result[i],v[1],v[2]);
			#nn[i]:=One(APR)*Value(nn[i],v[1],v[2]);
		od;
	od;
	return result;
end;