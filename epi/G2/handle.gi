
#
# Read("~/Workspace/groupsSB/epi/G2/handle.gi");
#

type:="G";
rank:=2;
nr_pos_roots:=6;

Read("~/Workspace/groupsSB/epi/group.gi");


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
