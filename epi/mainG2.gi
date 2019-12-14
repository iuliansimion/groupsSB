
#
# Read("~/Workspace/epimorphic/mainG2.gi");
#

type:="G";
rank:=2;
nr_pos_roots:=6;



ZZ:=Integers;
avarnames:=List([1..10],i->Concatenation("P_{",String(i),"}(a)"));
bvarnames:=List([1..10],i->Concatenation("P_{",String(i),"}(b)"));
cvarnames:=List([1..10],i->Concatenation("P_{",String(i),"}(c)"));
varnames:=Concatenation(avarnames,bvarnames,cvarnames);
APR:=PolynomialRing(ZZ,varnames);
vars:=IndeterminatesOfPolynomialRing(APR);

sla:=SimpleLieAlgebraTypeA_G(type,rank,APR);

cb:=CanonicalBasis(sla);

e:=cb[1];

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
	result:=tmp^0;
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

 #[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

handleUaUbUc:=function()
	local rels, vals,i,v;
	# UaUb=Uc = U(a1+b1,a2+b2,a3+b3+a2*b1)
	rels:=Set(Concatenation(Ua*Ub-Uc));
	
	vals:=[];
	Append(vals,[[[vars[22]],[vars[2]+vars[12]]]]);
	Append(vals,[[[vars[21]],[vars[1]+vars[11]]]]);
	Append(vals,[[[vars[23]],[vars[3]+vars[13]+vars[2]*vars[11]]]]);
#a_{2}^2*b_{1}+2*a_{2}*b_{1}*b_{2}-2*a_{3}*b_{2}-a_{4}-b_{4}
#-a_{2}*b_{1}^2-2*a_{3}*b_{1}-a_{4}-b_{4}+c_{4},
	Append(vals,[[[vars[24]],[vars[2]*vars[11]^2+2*vars[3]*vars[11]+vars[4]+vars[14]]]]);
#a_{2}*b_{1}^3+3*a_{3}*b_{1}^2+3*a_{4}*b_{1}+a_{5}+b_{5}-c_{5}
	Append(vals,[[[vars[25]],[vars[2]*vars[11]^3+3*vars[3]*vars[11]^2+3*vars[4]*vars[11]+vars[5]+vars[15]]]]);
#a_{2}^2*b_{1}^3
#-a_{2}*b_{1}^3*b_{2}
#+3*a_{2}*a_{3}*b_{1}^2
#+3*a_{2}*b_{1}^2*b_{3}
#-3*a_{3}*b_{1}^2*b_{2}

#+3*a_{3}^2*b_{1}
#+6*a_{3}*b_{1}*b_{3}
#-3*a_{4}*b_{1}*b_{2}
#+3*a_{4}*b_{3}
#-a_{5}*b_{2}
#-a_{6}-b_{6}+c_{6}
	Append(vals,[[[vars[26]],[
-vars[2]^2*vars[11]^3
+vars[2]*vars[11]^3*vars[12]
-3*vars[2]*vars[3]*vars[11]^2
-3*vars[2]*vars[11]^2*vars[13]
+3*vars[3]*vars[11]^2*vars[12]
-3*vars[3]^2*vars[11]
-6*vars[3]*vars[11]*vars[13]
+3*vars[4]*vars[11]*vars[12]
-3*vars[4]*vars[13]
+vars[5]*vars[12]
+vars[6]+vars[16]]]]);


	for i in [1..Length(rels)] do
		for v in vals do
			Print(rels[i],"\n");
			rels[i]:=One(APR)*Value(rels[i],v[1],v[2]);
		od;
	od;
	return [vals,rels];
end;

handleUaUaUc:=function()
	local rels, vals,i,v;
	# Ua^2=Uc = U(2*a1,2*a2,2*a3+a1*a2)
	rels:=Set(Concatenation(Ua*Ua-Uc));
	
	vals:=[];
	Append(vals,[[[vars[22]],[vars[2]+vars[2]]]]);
	Append(vals,[[[vars[21]],[vars[1]+vars[1]]]]);
	Append(vals,[[[vars[23]],[vars[3]+vars[3]+vars[2]*vars[1]]]]);


	for i in [1..Length(rels)] do
		for v in vals do
			Print(rels[i],"\n");
			rels[i]:=One(APR)*Value(rels[i],v[1],v[2]);
		od;
	od;
end;

#sList(rels,r->Value(r,[vars[22]],[vars[2]+vars[12]]);


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


#
# c_2=0 c_3=0
#
#
#
handle23:=function()
	ua:=root_group(1,1);
	uabb:=root_group(4,1);
	u:=ua*uabb;
	id_mat:=u^0;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
	Append(vals,[[[xvars[7]],[Zero(APR)]]]);
	Append(vals,[[[xvars[8]],[Zero(APR)]]]);
	Append(vals,[[[xvars[5]],[Zero(APR)]]]);
	Append(vals,[[[xvars[10]],[Zero(APR)]]]); # for p>2
	Append(vals,[[[xvars[9]],[Zero(APR)]]]); # for p>2
	Append(vals,[[[xvars[2]],[-xvars[6]]]]);
	# ###
	# Append(vals,[[[xvars[1]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[3]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[4]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[8]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[9]],[Zero(APR)]]]);


	#
	# atentie la ordinea transpunerilor
	#
	rels:=nn*TransposedMat(u)-nn;
	rels:=evaluate_rels(rels,vals);
	nn:=evaluate_rels(nn,vals);

	return [u,nn,rels,vals];
end;
