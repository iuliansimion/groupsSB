#
# Read("~/Workspace/groupsSB/epi/A2/extTorusCentralizer.gi");
#

type:="A";
rank:=2;
nr_pos_roots:=3;

Read("~/Workspace/groupsSB/epi/group.gi");

highest_weight:=[1,1];# corresponding highest weight module is the Lie algebra

Read("~/Workspace/groupsSB/epi/ext.gi");


#Vb:=Basis(V);

get_weights_on_ext:=function(k)
    local d,Vb,on_module,weights_module;
    d:=Dimension(extA[k]);
    Vb:=extA_basis[k];
    #
    # action of elements in the chevalley basis cb0 on module
    #
    on_module:=List(cb0,b->List(Vb,v->Coefficients(Vb,b^v)));

    #
    # Wights of torus cb0[7],cb0[8] on module
    #
    weights_module:=[];

    for i in [1..d] do
        Add(weights_module,[on_module[7][i][i],on_module[8][i][i]]);
    od;
    
    return weights_module;
end;

get_0_weight_basis_on_ext:=function(k)
    local weights_module,pos;
    Vb:=extA_basis[k];
    weights_module:=get_weights_on_ext(k);
    pos:=Positions(weights_module,[0,0]);
    return Vb{pos};
end;

#Add(weights_module,[on_module[7][2][2],on_module[8][2][2]]);
#Add(weights_module,[on_module[7][3][3],on_module[8][3][3]]);

#Add(weights_module,[on_module[7][4][4],on_module[8][4][4]]);
#Add(weights_module,[on_module[7][5][5],on_module[8][5][5]]);

#Add(weights_module,[on_module[7][6][6],on_module[8][6][6]]);
#Add(weights_module,[on_module[7][7][7],on_module[8][7][7]]);
#Add(weights_module,[on_module[7][8][8],on_module[8][8][8]]);

#
# action of elements in the chevalley basis cb0 on module
#
on_algebra:=List(cb0,b->List(cb0,v->Coefficients(cb0,v*b)));

#
# Let a=cb0[1] and b=cb0[2]
#

#
# Vb[1] corrsponds to a+b (heighest root) ---> cb0[3]
#

#
# Vb[2] and Vb[3] corrspond to a and b
#

#
# In the internal notation
#   y1 = cb0[4] corrseponding to -a corresponding to cb0[7]
#   y2 = cb0[5] corresponding to -b corresponding to cb0[8]
#   y3 = cb0[6] corresponding to -a-b
