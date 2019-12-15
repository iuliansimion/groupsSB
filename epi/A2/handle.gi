
#
# Read("~/Workspace/groupsSB/epi/A2/handle.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

Read("~/Workspace/groupsSB/epi/group.gi");



#sList(rels,r->Value(r,[vars[22]],[vars[2]+vars[12]]);

#
# c_1, c_2 not 0 => c_3=c_1/2
#
handleReg:=function()
	# a test that the coeffs are right
	#a:=root_group(1,xvars[1]);
	#b:=root_group(2,xvars[2]);
	#c:=root_group(3,xvars[3]);
	#test1:=evaluate_U(a*b*c,[[[xvars[1]],[One(APR)]],[[xvars[2]],[One(APR)]],[[xvars[3]],[One(APR)/2]]]);
	#test2:=evaluate_U(a*b*c,[[[xvars[1]],[One(APR)*2]],[[xvars[2]],[One(APR)*2]],[[xvars[3]],[One(APR)*2]]]);
	#test2-test1^2;

	local ua,ub,uab,u,e,nn,vals,rels,i,v;
	#ua:=root_group(1,2);
	#ub:=root_group(2,2);
	#uab:=root_group(3,2);
	ua:=root_group(1,1);
	ub:=root_group(2,1);
	uab:=root_group(3,1/2);
	u:=ua*ub*uab;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
	Append(vals,[[[xvars[6]],[Zero(APR)]]]);
	Append(vals,[[[xvars[4]],[Zero(APR)]]]);
	Append(vals,[[[xvars[5]],[Zero(APR)]]]);
	Append(vals,[[[xvars[7]],[2*xvars[8]]]]);
	Append(vals,[[[xvars[8]],[Zero(APR)]]]);
	Append(vals,[[[xvars[2]],[xvars[1]]]]);
	##
	#Append(vals,[[[xvars[1]],[Zero(APR)]]]);
	#Append(vals,[[[xvars[3]],[Zero(APR)]]]);
	#Append(vals,[[[xvars[7]],[Zero(APR)]]]);

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
	return [u,nn,rels];
end;

#
# c_2=0
#
handle2:=function()
	local ua,ub,uab,u,e,nn,vals,rels,i,v;
	ua:=root_group(1,1);
	uab:=root_group(3,1);
	u:=ua*uab;
	e:=u-id_mat;

	nn:=List([1..2*nr_pos_roots+rank],i->xvars[i]);

	vals:=[];
	Append(vals,[[[xvars[4]],[Zero(APR)]]]);
	Append(vals,[[[xvars[6]],[Zero(APR)]]]);
	Append(vals,[[[xvars[2]],[-xvars[7]-xvars[8]]]]);
	Append(vals,[[[xvars[5]],[-2*xvars[7]+xvars[8]]]]);
	##
	Append(vals,[[[xvars[1]],[Zero(APR)]]]);
	Append(vals,[[[xvars[3]],[Zero(APR)]]]);
	Append(vals,[[[xvars[7]],[Zero(APR)]]]);

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
	return [u,nn,rels];
end;


