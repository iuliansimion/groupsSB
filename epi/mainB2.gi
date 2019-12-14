
#
# Read("~/Workspace/epimorphic/mainB2.gi");
#

type:="B";
rank:=2;
nr_pos_roots:=4;



ZZ:=Integers;
avarnames:=List([1..10],i->Concatenation("a_{",String(i),"}"));
bvarnames:=List([1..10],i->Concatenation("b_{",String(i),"}"));
cvarnames:=List([1..10],i->Concatenation("c_{",String(i),"}"));
xvarnames:=List([1..10],i->Concatenation("x_{",String(i),"}"));
varnames:=Concatenation(avarnames,bvarnames,cvarnames,xvarnames);
APR:=PolynomialRing(ZZ,varnames);
vars:=IndeterminatesOfPolynomialRing(APR);
xvars:=vars{[31..40]};

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

handleUaUbUc:=function()
	local rels, vals,i,v;
	# UaUb=Uc = U(a1+b1,a2+b2,a3+b3+a2*b1)
	rels:=Set(Concatenation(Ua*Ub-Uc));
	
	vals:=[];
	Append(vals,[[[vars[22]],[vars[2]+vars[12]]]]);
	Append(vals,[[[vars[21]],[vars[1]+vars[11]]]]);
	Append(vals,[[[vars[23]],[vars[3]+vars[13]-vars[2]*vars[11]]]]);
#a_{2}^2*b_{1}+2*a_{2}*b_{1}*b_{2}-2*a_{3}*b_{2}-a_{4}-b_{4}
	Append(vals,[[[vars[24]],[-vars[2]^2*vars[11]-2*vars[2]*vars[11]*vars[12]+2*vars[3]*vars[12]+vars[4]+vars[14]]]]);


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

#sList(rels,r->Value(r,[vars[22]],[vars[2]+vars[12]]);

#
# c_1, c_2 not 0 => c_3=c_1/2
#
# a test that the coeffs are right
#a:=root_group(1,xvars[1]);
#b:=root_group(2,xvars[2]);
#c:=root_group(3,xvars[3]);
#d:=root_group(4,xvars[4]);
#test1:=evaluate_U(a*b*c*d,[[[xvars[1]],[One(APR)]],[[xvars[2]],[One(APR)]],[[xvars[3]],[-One(APR)/2]],[[xvars[4]],[-2*One(APR)/3]]]);
#test2:=evaluate_U(a*b*c*d,[[[xvars[1]],[One(APR)*2]],[[xvars[2]],[One(APR)*2]],[[xvars[3]],[-One(APR)*2]],[[xvars[4]],[-16*One(APR)/3]]]);
#Set(Concatenation(test2-test1^2));
handleReg:=function()
	local ua,ub,uab,u,e,nn,vals,rels,i,v;
	ua:=root_group(1,1);
	ub:=root_group(2,1);
	uab:=root_group(3,-1/2);
	uabb:=root_group(4,-2/3);
	u:=ua*ub*uab*uabb;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
 	Append(vals,[[[xvars[6]],[Zero(APR)]]]);
# 	Append(vals,[[[xvars[4]],[Zero(APR)]]]);
  	Append(vals,[[[xvars[5]],[Zero(APR)]]]);
  	Append(vals,[[[xvars[7]],[Zero(APR)]]]);
  	Append(vals,[[[xvars[8]],[Zero(APR)]]]);
  	Append(vals,[[[xvars[9]],[2*xvars[10]]]]);
 	Append(vals,[[[xvars[10]],[Zero(APR)]]]);
  	Append(vals,[[[xvars[2]],[xvars[1]]]]);
  	Append(vals,[[[xvars[3]],[Zero(APR)]]]);
#  	Append(vals,[[[xvars[3]],[-2*xvars[1]/3]]]);
#  	##
	#Append(vals,[[[xvars[1]],[Zero(APR)]]]);
	#Append(vals,[[[xvars[3]],[Zero(APR)]]]);
	#Append(vals,[[[xvars[7]],[Zero(APR)]]]);

	#
	# atentie la ordinea transpunerilor
	#
	rels:=nn*TransposedMat(u)-nn;
	rels:=evaluate_rels(rels,vals);
	nn:=evaluate_rels(nn,vals);
	
	return [u,nn,rels,vals];
end;

#
# Done
#
handleRegModule:=function()
	local ua,ub,uab,u,e,nn,vals,rels,i,v;
	ua:=root_group(1,1);
	ub:=root_group(2,1);
	uab:=root_group(3,-1/2);
	uabb:=root_group(4,-2/3);
	u:=ua*ub*uab*uabb;
	
	nn:=One(APR)*[0,0,0,0,xvars[1],xvars[2],xvars[3],0,0,0];
	#
	# atentie la ordinea transpunerilor
	#
	nn:=nn*TransposedMat(u);
	
	rels:=[];
	Append(rels,[nn[1]-nn[2]]);
	# Append(rels,[nn[3]-2*nn[1]/3]);
	
	vals:=[];
	Append(vals,[[[xvars[3]],[-6*xvars[1]/5+8*xvars[2]/5]]]); # because I want the projection on U3 to be trivial
	Append(vals,[[[xvars[2]],[3*xvars[1]/4]]]); # because I want the projections on U1 and U2 to be equal
#	Append(vals,[[[xvars[3]],[13*xvars[1]/2-2*xvars[2]]]]);
#	Append(vals,[[[xvars[2]],[203*3*xvars[1]/(4*85)]]]);

	rels:=evaluate_rels(rels,vals);
	nn:=evaluate_rels(nn,vals);
	
	return [u,nn,rels,vals];
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

#
# c_2=0 c_4=0
#
handle24:=function()
	ua:=root_group(1,1);
	uab:=root_group(3,1);
	u:=ua*uab;
	id_mat:=u^0;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
	Append(vals,[[[xvars[7]],[Zero(APR)]]]);
	Append(vals,[[[xvars[8]],[Zero(APR)]]]);
	Append(vals,[[[xvars[5]],[Zero(APR)]]]);
	Append(vals,[[[xvars[2]],[Zero(APR)]]]); # for p>2
	Append(vals,[[[xvars[9]],[Zero(APR)]]]);
	Append(vals,[[[xvars[10]],[-xvars[6]]]]); # for p>2
	# # ###
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


#
# c_2=0
#
handle2:=function()
	ua:=root_group(1,1);
	uab:=root_group(3,1);
	uabb:=root_group(4,1);
	u:=ua*uab*uabb;
	id_mat:=u^0;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
	Append(vals,[[[xvars[5]],[xvars[8]]]]);
	Append(vals,[[[xvars[7]],[-xvars[8]]]]);
	Append(vals,[[[xvars[2]],[-xvars[6]+xvars[9]]]]);
	Append(vals,[[[xvars[6]],[xvars[9]-xvars[10]]]]); # only in char 2
	###
	# Append(vals,[[[xvars[1]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[3]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[4]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[8]],[Zero(APR)]]]);
	# Append(vals,[[[xvars[9]],[Zero(APR)]]]);


	#
	# atentie la ordinea transpunerilor
	#
	rels:=nn*TransposedMat(u)-nn;
	for i in [1..Length(rels)] do
		for v in vals do
			Print(rels[i],"\n");
			rels[i]:=One(APR)*Value(rels[i],v[1],v[2]);
			nn[i]:=One(APR)*Value(nn[i],v[1],v[2]);
		od;
	od;

	return [u,nn,rels,vals];
end;


#
# regular
# c_1, c_2 not, c_3 not 0
#
handle1:=function()

end;

