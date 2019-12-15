
#
# Read("~/Workspace/groupsSB/epi/A2/generic.UU.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

Read("~/Workspace/groupsSB/epi/group.gi");

handleUaUbUc:=function()
	local rels, vals,i,v;
	# UaUb=Uc = U(a1+b1,a2+b2,a3+b3+a2*b1)
	rels:=Set(Concatenation(Ua*Ub-Uc));
	
	vals:=[];
	Append(vals,[[[vars[22]],[vars[2]+vars[12]]]]);
	Append(vals,[[[vars[21]],[vars[1]+vars[11]]]]);
	Append(vals,[[[vars[23]],[vars[3]+vars[13]+vars[2]*vars[11]]]]);


	for i in [1..Length(rels)] do
		for v in vals do
			Print(rels[i],"\n");
			rels[i]:=One(APR)*Value(rels[i],v[1],v[2]);
		od;
	od;
	return vals;
end;
