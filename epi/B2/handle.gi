
#
# Read("~/Workspace/groupsSB/epi/B2/handle.gi");
#

type:="B";
rank:=2;
nr_pos_roots:=4;

Read("~/Workspace/groupsSB/epi/group.gi");


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
	local ua,ub,uab,uabb,u,e,nn,vals,rels,i,v;
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
	local ua,ub,uab,uabb,u,e,nn,vals,rels,i,v;
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
	local ua,ub,uab,uabb,u,e,nn,vals,rels,i,v;
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
	local ua,ub,uab,uabb,u,e,nn,vals,rels,i,v;
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
	local ua,ub,uab,uabb,u,e,nn,vals,rels,i,v;
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

