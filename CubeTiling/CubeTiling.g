# we search for a counter example to Keller conjecture:
# every cube tiling of R^n has two cubes with two facets
# in common
# we use the formalism of Lagarias Shor using cliques.






SymmetryGroup:=function(n)
  local VertexSet, ListGen, ListMat, TheMat, ListV, i, eV, iP, g, j, GrpPerm, GrpMat, phi;
  VertexSet:=BuildSet(n, [0,1,2,3]);
  Print("VertexSet Generated\n");

  ListGen:=[];
  ListMat:=[];
  #
  ListV:=[];
  for i in [1..Length(VertexSet)]
  do
    eV:=ShallowCopy(VertexSet[i]);
    eV[1]:=(eV[1]+1) mod 4;
    Add(ListV, Position(VertexSet, eV));
  od;
  TheMat:=IdentityMat(n+1);
  TheMat[n+1][1]:=1/4;
  Add(ListGen, PermList(ListV));
  Add(ListMat, TheMat);
  #
  ListV:=[];
  for i in [1..Length(VertexSet)]
  do
    eV:=ShallowCopy(VertexSet[i]);
    eV[1]:=3-eV[1];
    Add(ListV, Position(VertexSet, eV));
  od;
  TheMat:=IdentityMat(n+1);
  TheMat[1][1]:=-1;
  TheMat[n+1][1]:=1;
  Add(ListGen, PermList(ListV));
  Add(ListMat, TheMat);
  #
  for iP in [1..n-1]
  do
    g:=(iP, iP+1);
    ListV:=[];
    for i in [1..Length(VertexSet)]
    do
      ListV[i]:=Position(VertexSet, Permuted(VertexSet[i], g));
    od;
    TheMat:=NullMat(n+1, n+1);
    TheMat[n+1][n+1]:=1;
    for i in [1..n]
    do
      TheMat[i][OnPoints(i, g)]:=1;
    od;
    Add(ListGen, PermList(ListV));
    Add(ListMat, TheMat);
  od;
  GrpPerm:=Group(ListGen);
  GrpMat:=Group(ListMat);
  phi:=GroupHomomorphismByImages(GrpPerm, GrpMat, ListGen, ListMat);
  Print("Symmetry group created\n");
  return rec(VertexSet:=VertexSet, SymmetryGroup:=GrpPerm, phi:=phi, GrpMat:=GrpMat, ListGen:=ListGen, ListMat:=ListMat);
end;



IsAdjacentStarGraph:=function(m1, m2)
  local nbDiff, i, testM;
  nbDiff:=0;
  testM:=false;
  for i in [1..Length(m1)]
  do
    if m1[i]<>m2[i] then
      nbDiff:=nbDiff+1;
    fi;
    if AbsInt(m1[i]-m2[i])=2 then
      testM:=true;
    fi;
  od;
  if testM=true and nbDiff>=2 then
    return true;
  else
    return false;
  fi;
end;




IsAdjacent:=function(m1, m2)
  local i;
  for i in [1..Length(m1)]
  do
    if AbsInt(m1[i]-m2[i])=2 then
      return true;
    fi;
  od;
  return false;
end;



BuildUnextendible:=function(EXMP, p)
  local NewClique, eV, fV, V, VSET;
  NewClique:=[];
  VSET:=BuildSet(p, [0,2]);
  for eV in EXMP
  do
    for fV in VSET
    do
      V:=ShallowCopy(eV);
      Append(V, fV);
      Add(NewClique, V);
    od;
  od;
  return NewClique;
end;



FuncListMatching:=function(LS)
  local n, iCol, ListMatch, eV, val;
  n:=Length(LS[1]);
  for iCol in [1..n]
  do
    ListMatch:=[0,0,0,0];
    for eV in LS
    do
      val:=eV[iCol];
      ListMatch[val+1]:=ListMatch[val+1]+1;
    od;
  od;
  Print("0:", ListMatch[1],"  1:", ListMatch[2],"  2:", ListMatch[3],"  3:", ListMatch[4], "\n");
end;


NzFunctionBis:=function(Clique, z)
  local n, Nz, eVect, eDiff, test, iCol;
  n:=Length(Clique[1]);
  Nz:=0;
  for eVect in Clique
  do
    eDiff:=List(eVect-z, x-> x mod 4);
    test:=1;
    for iCol in [1..n]
    do
      if eDiff[iCol]=3 then
        test:=0;
      fi;
    od;
    if test=1 then
      Nz:=Nz+1;
    fi;
  od;
  return Nz;
end;






Is2ZnPeriodic:=function(eClique)
  local eSet, ListVal, eVect, z, eVal, n;
  n:=Length(eClique[1]);
  for eSet in BuildSet(n, [0,1])
  do
    ListVal:=[];
    for eVect in BuildSet(n, [0,2])
    do
      z:=eSet+eVect;
      eVal:=NzFunctionBis(eClique, z);
      AddSet(ListVal, eVal);
      if Length(ListVal) > 1 then
        return false;
      fi;
    od;
  od;
  return true;
end;


PossibilityOfExtension:=function(Clique, VSET)
  local n, FuncTest, eV, ListMatch;
  n:=Length(Clique[1]);
  FuncTest:=function(eV)
    local eC;
    for eC in Clique
    do
      if IsAdjacent(eV, eC)=false then
        return false;
      fi;
    od;
    return true;
  end;
  ListMatch:=[];
  for eV in BuildSet(n, VSET)
  do
    if FuncTest(eV)=true then
      Add(ListMatch, eV);
    fi;
  od;
  return ListMatch;
end;


IsCubePackingDetermining:=function(Clique, VSET)
  local n, ListPoss, nbPoss, iV, jV;
  n:=Length(Clique[1]);
  ListPoss:=PossibilityOfExtension(Clique, VSET);
  nbPoss:=Length(ListPoss);
  for iV in [1..nbPoss-1]
  do
    for jV in [iV+1..nbPoss]
    do
      if IsAdjacent(ListPoss[iV], ListPoss[jV])=false then
        return false;
      fi;
    od;
  od;
  return true;
end;




IsExtendible:=function(Clique)
  local ListExt;
  ListExt:=PossibilityOfExtension(Clique, [0,1,2,3]);
  if Length(ListExt)=0 then
    return false;
  else
    return true;
  fi;
end;



SearchPossibilityExtension:=function(Clique, VSET)
  local n, eV, FuncTest;
  n:=Length(Clique[1]);
  FuncTest:=function(eV)
    local eC;
    for eC in Clique
    do
      if IsAdjacent(eV, eC)=false then
        return false;
      fi;
    od;
    return true;
  end;
  for eV in BuildSet(n, VSET)
  do
    if FuncTest(eV)=true then
      Print("eV=", eV, "\n");
      return true;
    fi;
  od;
  return false;
end;





RandomElement:=function(n, VSET)
  local eV, i;
  eV:=[];
  for i in [1..n]
  do
    Add(eV, Random(VSET));
  od;
  return eV;
end;


TestTilingNess:=function(Clique)
  local n, i, j;
  n:=Length(Clique[1]);
  if Length(Clique)<>2^n then
    return false;
  fi;
  for i in [1..Length(Clique)-1]
  do
    for j in [i+1..Length(Clique)]
    do
      if IsAdjacent(Clique[i], Clique[j])=false then
        return false;
      fi;
    od;
  od;
  return true;
end;



MaxSearch:=2000;
RandomSearch:=function(n, AskedSize)
  local Clique, iter, nbIter, SUC, RND, iPos, test;
  Clique:=[ShallowCopy(ListWithIdenticalEntries(n, 0))];
  for iter in [1..AskedSize]
  do
    nbIter:=1;
    SUC:=false;
    while(nbIter<MaxSearch)
    do
      RND:=RandomElement(n, [0,1,2,3]);
      iPos:=1;
      test:=true;
      while(iPos<=Length(Clique))
      do
        if IsAdjacentStarGraph(RND, Clique[iPos])=false then
          test:=false;
          break;
        fi;
        iPos:=iPos+1;
      od;
      if test=true then
        Add(Clique, RND);
        SUC:=true;
        break;
      fi;
      nbIter:=nbIter+1;
    od;
#    Print("We have ", Length(Clique), " elements\n");
    if SUC=false then
      return Clique;
    fi;
  od;
  return Clique;
end;








#
#
# if VSET=[0,1,2,3] then we do random tiling
# while if VSET=[0,1,2] we do random packing.
RandomPackingTiling:=function(n, VSET)
  local Clique, iter, SUC, RND, iPos, test;
  Clique:=[ShallowCopy(ListWithIdenticalEntries(n, 0))];
  while(true)
  do
    if SearchPossibilityExtension(Clique, VSET)=false then
      break;
    fi;
    while(true)
    do
      RND:=RandomElement(n, VSET);
      iPos:=1;
      test:=true;
      while(iPos<=Length(Clique))
      do
        if IsAdjacent(RND, Clique[iPos])=false then
          test:=false;
          break;
        fi;
        iPos:=iPos+1;
      od;
      if test=true then
        Add(Clique, RND);
        break;
      fi;
    od;
#    Print("We have ", Length(Clique), " elements\n");
  od;
  return Clique;
end;





FindNumberTiles:=function(Clique, ListC)
  local nbr, eTile, Match, ePos;
  nbr:=0;
  for eTile in Clique
  do
    Match:=true;
    for ePos in ListC
    do
      if not(eTile[ePos] in [0,1]) then
        Match:=false;
      fi;
    od;
    if Match=true then
      nbr:=nbr+1;
    fi;
  od;
  return nbr;
end;


#
# useful only for partial cliques
ReducedSymbol:=function(Clique, z)
  local n, u, FuncCompute, ListMatch, val;
  n:=Length(Clique[1]);
  FuncCompute:=function(eV)
    local h;
    for h in BuildSet(n, [0,1])
    do
      if Position(Clique, (h+eV) mod 4)<>fail then
        return 1;
      fi;
    od;
    return 0;
  end;
  ListMatch:=[];
  for u in BuildSet(n, [0,2])
  do
    val:=FuncCompute(z+u);
    if val=1 then
      Add(ListMatch, z+u);
    fi;
  od;
  return ListMatch;
end;


IsAffineSet:=function(ListVect)
  local NewSet, eV, u, v;
  NewSet:=[];
  for eV in ListVect
  do
    Add(NewSet, (eV-ListVect[1]) mod 4);
  od;
  for u in NewSet
  do
    for v in NewSet
    do
      if Position(NewSet, (u+v) mod 4)=fail then
        return false;
      fi;
    od;
  od;
  return true;
end;



AssociatedWithAffineSet:=function(Clique)
  local n, ListMatch, h, ERL, test;
  n:=Length(Clique[1]);
  ListMatch:=[];
  for h in BuildSet(n, [0,1])
  do
    ERL:=ReducedSymbol(Clique, h);
    test:=IsAffineSet(ERL);
    if test=true then
      Add(ListMatch, h);
    fi;
  od;
  return ListMatch;
end;




RandomEquivalenceOfClique:=function(Clique)
  local n, ePerm, eSign, eAdd, i, eTile, fTile, gTile, NewClique;
  n:=Length(Clique[1]);
  ePerm:=Random(SymmetricGroup(n));
  eSign:=[];
  for i in [1..n]
  do
    eSign[i]:=Random([0,1]);
  od;
  eAdd:=[];
  for i in [1..n]
  do
    eAdd[i]:=Random([0,1,2,3]);
  od;
  NewClique:=[];
  for eTile in Clique
  do
    fTile:=Permuted(eTile, ePerm);
    gTile:=[];
    for i in [1..n]
    do
      if eSign[i]=0 then
        gTile[i]:=(fTile[i]+eAdd[i]) mod 4;
      else
        gTile[i]:=((3-fTile[i])+eAdd[i]) mod 4;
      fi;
    od;
    Add(NewClique, gTile);
  od;  
  return NewClique;
end;



TwiceSymGrp:=function(n)
  local ListGen, i;
  ListGen:=[];
  for i in [1..n]
  do
    Add(ListGen, (2*i-1, 2*i));
  od;
  for i in [1..n-1]
  do
    Add(ListGen, (2*i-1, 2*i+1)(2*i, 2*i+2));
  od;
  return Group(ListGen);
end;

TwiceRepresentation:=function(eV)
  local NewV, i, h;
  NewV:=[];
  for i in [1..Length(eV)]
  do
    Append(NewV, [eV[i],(4-eV[i]) mod 4]);
  od;
  return NewV;
end;

BackFromTwiceRepresentation:=function(eV)
  local n, i, NewV;
  n:=Length(eV)/2;
  NewV:=[];
  for i in [1..n]
  do
    NewV[i]:=eV[2*i-1];
  od;
  return NewV;
end;


CreateVectorSet:=function(Clique, idx)
  local ListVect, iTile, NewV;
  ListVect:=[];
  for iTile in [1..Length(Clique)]
  do
    if iTile<>idx then
      NewV:=TwiceRepresentation(Clique[iTile]-Clique[idx]);
      AddSet(ListVect, NewV);
    fi;
  od;
  return ListVect;
end;


ActionOnVectSet:=function(ListVect, g)
  local NewSet, eV;
  NewSet:=[];
  for eV in ListVect
  do
    AddSet(NewSet, Permuted(eV, g));
  od;
  return NewSet;
end;



IsEquivalentClique:=function(Clique1, Clique2)
  local n, SYM, idx, ListVect1, ListVect2, g;
  n:=Length(Clique1[1]);
  SYM:=TwiceSymGrp(n);
  ListVect1:=CreateVectorSet(Clique1, 1);
  for idx in [1..Length(Clique2)]
  do
    ListVect2:=CreateVectorSet(Clique2, idx);
    g:=RepresentativeAction(SYM, ListVect1, ListVect2, ActionOnVectSet);
    if g<>fail then
      return true;
    fi;
  od;
  return false;
end;




StabilizerOfPointInClique:=function(Clique, idx)
  local n, ListGen, ListV, iTile, SYM, ListVect, Stab, eGen, NewPerm, fV;
  n:=Length(Clique[1]);
  ListV:=[];
  for iTile in [1..Length(Clique)]
  do
    Add(ListV, TwiceRepresentation(Clique[iTile]-Clique[idx]));
  od;
  SYM:=TwiceSymGrp(n);
  ListVect:=CreateVectorSet(Clique, 1);
  Stab:=Stabilizer(SYM, ListVect, ActionOnVectSet);
  ListGen:=[];
  for eGen in GeneratorsOfGroup(Stab)
  do
    NewPerm:=[];
    for iTile in [1..Length(Clique)]
    do
      fV:=Permuted(ListV[iTile], eGen);
      NewPerm[iTile]:=Position(ListV, fV);
    od;
    Add(ListGen, PermList(NewPerm));
  od;
  return ListGen;
end;


InvariantOfClique:=function(Clique)
  local n, ListInv, i, j, V, iCol;
  n:=Length(Clique[1]);
  ListInv:=[];
  for i in [1..Length(Clique)-1]
  do
    for j in [i+1..Length(Clique)]
    do
      V:=(Clique[i]-Clique[j]) mod 4;
      for iCol in [1..n]
      do
        if V[iCol]=3 then
          V[iCol]:=1;
        fi;
      od;
      Sort(V);
      Add(ListInv, ShallowCopy(V));
    od;
  od;
  return Collected(ListInv);
end;







SymmetryOfTiling:=function(Clique)
  local GRA, i, j, n, SymGrp, MatDist, ListCase, CHR, eCol, DE, eCase;
  n:=Length(Clique[1]);
  GRA:=NullGraph(Group(()), Length(Clique));
  for i in [1..Length(Clique)-1]
  do
    for j in [i+1..Length(Clique)]
    do
      if IsAdjacentStarGraph(Clique[i], Clique[j])=true then
        AddEdgeOrbit(GRA, [i,j]);
        AddEdgeOrbit(GRA, [j,i]);
      fi;
    od;
  od;
  SymGrp:=AutGroupGraph(GRA);
  MatDist:=NullMat(Length(Clique), Length(Clique));
  ListCase:=[];
  for i in [1..Length(Clique)-1]
  do
    for j in [i+1..Length(Clique)]
    do
      CHR:=[];
      for eCol in [1..n]
      do
        DE:=Clique[i][eCol]-Clique[j][eCol];
        if DE>=3 then 
          DE:=DE-4;
        fi;
        if DE<=-3 then 
          DE:=DE+4;
        fi;
        Add(CHR, AbsInt(DE));
      od;
      Sort(CHR);
      MatDist[i][j]:=ShallowCopy(CHR);
      MatDist[j][i]:=ShallowCopy(CHR);
      AddSet(ListCase, ShallowCopy(CHR));
    od;
  od;
  for eCase in ListCase
  do
    GRA:=NullGraph(Group(()), Length(Clique));
    for i in [1..Length(Clique)-1]
    do
      for j in [i+1..Length(Clique)]
      do
        if MatDist[i][j]=eCase then
          AddEdgeOrbit(GRA, [i,j]);
          AddEdgeOrbit(GRA, [j,i]);
        fi;
      od;
    od;
    SymGrp:=Intersection(SymGrp, AutGroupGraph(GRA));
  od;
  return SymGrp;
end;



RegularTiling:=function(n)
  return BuildSet(n, [0,2]);
end;

MergeTiling:=function(Clique1, Clique2)
  local NewClique, eV, h;
  NewClique:=[];
  for eV in Clique1
  do
    h:=ShallowCopy(eV);
    Add(h, 0);
    Add(NewClique, h);
  od;
  for eV in Clique2
  do
    h:=ShallowCopy(eV);
    Add(h, 2);
    Add(NewClique, h);
  od;
  return NewClique;
end;





FacetMovements:=function(Clique)
  local n, nbElt, ListNewClique, i, j, eTile, fTile, ListDiff, iCol, eCol, NewClique, u;
  n:=Length(Clique[1]);
  nbElt:=2^n;
  ListNewClique:=[];
  for i in [1..nbElt-1]
  do
    for j in [i+1..nbElt]
    do
      eTile:=Clique[i];
      fTile:=Clique[j];
      ListDiff:=[];
      for iCol in [1..n]
      do
        if eTile[iCol]<>fTile[iCol] then
          Add(ListDiff, iCol);
        fi;
      od;
      if Length(ListDiff)=1 then
        eCol:=ListDiff[1];
        NewClique:=[];
        for u in [1..nbElt]
        do
          if u<>i and u<>j then
            Add(NewClique, Clique[u]);
          else
            eTile:=[];
            for iCol in [1..n]
            do
              if iCol<>eCol then
                Add(eTile, Clique[u][iCol]);
              else
                Add(eTile, (Clique[u][iCol]+1) mod 4);
              fi;
            od;
            Add(NewClique, eTile);
          fi;
        od;
        Add(ListNewClique, NewClique);
      fi;
    od;
  od;
  return ListNewClique;
end;





SearchColumnsOfColumns:=function(Clique)
  local n, FuncTestMatch, FuncTestDetermineColumn, eTile, fTile, ListOfMaximalColumn, eSet, ListCoord, MatchedSet, iSize, reply, i, j;
  n:=Length(Clique[1]);
  FuncTestMatch:=function(ListSet, eSet)
    local eV;
    for eV in ListSet
    do
      if IsSubset(eV, eSet)=true then
        return true;
      fi;
    od;
    return false;
  end;

  FuncTestDetermineColumn:=function(ListC, eTile)
    local eSet, fTile, i, ListTile;
    ListTile:=[];
    for eSet in Combinations(ListC)
    do
      fTile:=[];
      for i in [1..n]
      do
        if i in eSet then
          Add(fTile, (eTile[i]+2) mod 4);
        else
          Add(fTile, eTile[i]);
        fi;
      od;
      if Position(Clique, fTile)=fail then
        return false;
      else
        AddSet(ListTile, fTile);
      fi;
    od;
    return ListTile;
  end;
  ListOfMaximalColumn:=[];
  for eTile in Clique
  do
    ListCoord:=[];
    for i in [1..n]
    do
      fTile:=[];
      for j in [1..n]
      do
        if j<>i then
          Add(fTile, eTile[j]);
        else
          Add(fTile, (eTile[i]+2) mod 4);
        fi;
      od;
      if Position(Clique, fTile)<>fail then
        Add(ListCoord, i);
      fi;
    od;
    MatchedSet:=[];
    for iSize in Reversed([1..Length(ListCoord)])
    do
      for eSet in Combinations(ListCoord, iSize)
      do
        if FuncTestMatch(MatchedSet, eSet)=false then
          reply:=FuncTestDetermineColumn(eSet, eTile);
          if reply<>false then
            Add(MatchedSet, eSet);
            AddSet(ListOfMaximalColumn, reply);
          fi;
        fi;
      od;
    od;
  od;
  return ListOfMaximalColumn;
end;




SearchLayers:=function(Clique, level)
  local n, TestLayerness, ListLayer, eSet;
  n:=Length(Clique[1]);
  TestLayerness:=function(eSet)
    local ListProj, FuncInsert, eTile, TileProj, reply;
    ListProj:=[];
    FuncInsert:=function(eP)
      local eC;
      if Position(ListProj, eP)<>fail then
        return true;
      fi;
      for eC in ListProj
      do
        if IsAdjacent(eC, eP)=false then
          return false;
        fi;
      od;
      AddSet(ListProj, eP);
      return true;
    end;
    for eTile in Clique
    do
      TileProj:=eTile{eSet};
      reply:=FuncInsert(TileProj);
      if reply=false then
        return false;
      fi;      
    od;
    return true;
  end;

  ListLayer:=[];
  for eSet in Combinations([1..n], level)
  do
    if TestLayerness(eSet)=true then
      Add(ListLayer, eSet);
    fi;
  od;
  return ListLayer;
end;






NzFunction:=function(Clique, z, m)
  local n, ListXminimal, ListXmaximal, i, NSum, Prod, eV, u, v;
  n:=Length(Clique[1]);
  ListXminimal:=[];
  ListXmaximal:=[];
  for i in [1..n]
  do
    ListXminimal[i]:=z[i];
    ListXmaximal[i]:=z[i]+m;
  od;
  NSum:=0;
  for eV in Clique
  do
    Prod:=1;
    for i in [1..n]
    do
      u:=LowerInteger((ListXminimal[i]-eV[i])/4);
      v:=UpperInteger((ListXmaximal[i]-2-eV[i])/4);
      Prod:=Prod*(1+v-u);
    od;
    NSum:=NSum+Prod;
  od;
  return NSum;
end;




ComputeMoment:=function(Clique, m, k)
  local n, NSum, z, h;
  n:=Length(Clique[1]);
  NSum:=0;
  for z in BuildSet(n, [0,1,2,3])
  do
    h:=NzFunction(Clique, z, m);
#    Print("h=", h, "\n");
    NSum:=NSum+h^k;
  od;
  return NSum/(4^n);
end;





CTOrbitGroupFormalism:=function(nbE, GroupExt, ListOrbit)
  local Odisc, O, FuncInvariant, FuncInsert;
  Odisc:=[];
  O:=Orbits(GroupExt, Combinations([1..nbE], 1), OnSets);
  Append(Odisc, O);
  if Length(Odisc)<10 then
    O:=Orbits(GroupExt, Combinations([1..nbE], 2), OnSets);
    Append(Odisc, O);
  fi;
  if Length(Odisc)<10 then
    O:=Orbits(GroupExt, Combinations([1..nbE], 3), OnSets);
    Append(Odisc, O);
  fi;


  FuncInvariant:=function(ListInc)
    local eInv, eO, nb, eSet;
    eInv:=[];
    for eO in Odisc
    do
      nb:=0;
      for eSet in eO
      do
        if IsSubset(ListInc, eSet)=true then
          nb:=nb+1;
        fi;
      od;
      Add(eInv, nb);
    od;
    return eInv;
  end;


  FuncInsert:=function(Linc)
    local Stab, iOrb, TheInv, TheStab, iExt;
    TheInv:=FuncInvariant(Linc);
    for iOrb in [1..Length(ListOrbit)]
    do
      if TheInv=ListOrbit[iOrb].TheInv then
        if RepresentativeAction(GroupExt, Linc, ListOrbit[iOrb].Inc, OnSets)<>fail then
          return;
        fi;
      fi;
    od;
    Stab:=Stabilizer(GroupExt, Linc, OnSets);
    Add(ListOrbit, rec(Inc:=Linc, TheInv:=TheInv, Status:="NO", OrdStab:=Order(Stab)));
  end;
  return rec(FuncInsert:=FuncInsert);
end;


AddOneLevelOfSearch:=function(eCand, GRP, VertexSet, FuncAdj, FuncStor)
  local Cand, eOrb, Stab, eVert, AcceptedSet, FuncTest, fCand, ListCand, ListStor;
  FuncTest:=function(eVert)
    local fVert;
    for fVert in eCand
    do
      if FuncAdj(VertexSet[eVert], VertexSet[fVert])=false then
        return false;
      fi;
    od;
    return true;
  end;

  AcceptedSet:=[];
  for eVert in [1..Length(VertexSet)]
  do
    if Position(eCand, eVert)=fail then
      if FuncTest(eVert)=true then
        Add(AcceptedSet, eVert);
      fi;
    fi;
  od;
  Stab:=Stabilizer(GRP, eCand, OnSets);
  ListCand:=[];
  ListStor:=[];
  for eOrb in Orbits(Stab, AcceptedSet, OnPoints)
  do
    fCand:=Union(eCand, [eOrb[1]]);
    Add(ListCand, fCand);
    if FuncStor(VertexSet{fCand})=true then
      Add(ListStor, fCand);
    fi;
  od;
  return rec(ListCand:=ListCand, ListStor:=ListStor);
end;



OneStepEnumeration:=function(n, FuncAdj, ListInput, FuncStor)
  local VE, VertexSet, SymGrp, i, NewList, RPLcand, RPLstor, ListBlock, eCand, EVC, ListRed, U, eC, ListStor, ListStorRed;
  VE:=SymmetryGroup(n);
  VertexSet:=VE.VertexSet;
  SymGrp:=VE.SymmetryGroup;
  NewList:=[];
  ListBlock:=[];
  ListStor:=[];
  RPLcand:=CTOrbitGroupFormalism(Length(VertexSet), SymGrp, NewList);
  RPLstor:=CTOrbitGroupFormalism(Length(VertexSet), SymGrp, ListStor);
  for eCand in ListInput
  do
    EVC:=AddOneLevelOfSearch(eCand, SymGrp, VertexSet, FuncAdj, FuncStor);
    if Length(EVC.ListCand)=0 then
      Add(ListBlock, eCand);
    else
      for eC in EVC.ListCand
      do
        RPLcand.FuncInsert(eC);
      od;
    fi;
    for eC in EVC.ListStor
    do
      RPLstor.FuncInsert(eC);
    od;
  od;
  ListRed:=List(NewList, x->x.Inc);
  ListStorRed:=List(ListStor, x->x.Inc);
  return rec(ListRed:=ListRed, ListBlock:=ListBlock, ListStor:=ListStorRed);
end;




ExhaustiveEnumeration:=function(n, FuncAdj, FuncKill)
  local VE, VertexSet, SymGrp, i, Reord, LM, ListCand, eCand;
  VE:=SymmetryGroup(n);
  VertexSet:=VE.VertexSet;
  SymGrp:=VE.SymmetryGroup;
  ListCand:=[[1]];
  Print("We have one clique at step 1\n");
  for i in [2..2^n]
  do
    ListCand:=AddOneLevelOfSearch(ListCand, SymGrp, VertexSet, FuncAdj, FuncKill);
    Print("We have ", Length(ListCand), " cliques at step ", i, "\n");
  od;
  Reord:=[];
  for eCand in ListCand
  do
    LM:=[];
    for i in eCand
    do
      Add(LM, VertexSet[i]);
    od;
    Add(Reord, LM);
  od;
  return Reord;
end;



AdjacencyMethodTiling:=function(n)
  local VE, VertexSet, SymGrp, nbC, ListClique, ListDone, iCliq, fCliq, eCand, eV, eMin, fCand, FindUndone;
  VE:=SymmetryGroup(n);
  VertexSet:=VE.VertexSet;
  SymGrp:=VE.SymmetryGroup;
  ListClique:=[RegularTiling(n)];
  ListDone:=[0];
  while(true)
  do
    nbC:=Length(ListClique);
    for iCliq in [1..nbC]
    do
      if ListDone[iCliq]=0 then
        ListDone[iCliq]:=1;
        for fCliq in FacetMovements(ListClique[iCliq])
        do
          eCand:=[];
          for eV in fCliq
          do
            AddSet(eCand, Position(VertexSet, eV));
          od;
          eMin:=Minimum(Orbit(SymGrp, eCand, OnSets));
          fCand:=[];
          for eV in eMin
          do
            AddSet(fCand, VertexSet[eV]);
          od;
          if Position(ListClique, fCand)=fail then
            Add(ListClique, fCand);
            Add(ListDone, 0);
          fi;
        od;
      fi;
    od;
    Print("we have now ", Length(ListClique), " orbits of cliques\n");
    FindUndone:=false;
    for iCliq in [1..Length(ListDone)]
    do
      if ListDone[iCliq]=0 then
        FindUndone:=true;
      fi;
    od;
    if FindUndone=false then
      return ListClique;
    fi;
  od;
end;



MinimizeSecondMoment:=function(Clique)
  local n, TheClique, TheMoment, FindLower, fCliq, Mom;
  TheClique:=ShallowCopy(Clique);
  TheMoment:=ComputeMoment(TheClique, 4, 2);
  while(true)
  do
    FindLower:=false;
    Print("Beginning a new round\n");
    for fCliq in FacetMovements(TheClique)
    do
      Mom:=ComputeMoment(fCliq, 4, 2);
      if Mom<TheMoment then
        TheClique:=fCliq;
        TheMoment:=Mom;
        FindLower:=true;
        Print("Find a moment of ", TheMoment, "\n");
      fi;
    od;
    if FindLower=false then
      break;
    fi;
  od;
  return TheClique;
end;


DetectLocalizeHole:=function(Clique)
  local n, i, LS;
  n:=Length(Clique[1]);
  for i in [1..n]
  do
    LS:=Set(List(Clique, x->x[i]));
    if Length(LS)=2 then
      return true;
    fi;
  od;
  return false;
end;



ListVerticesUncovered:=function(Clique)
  local n, ListVert, IsCovered, eVert;
  n:=Length(Clique[1]);
  ListVert:=[];
  IsCovered:=function(eVect)
    local fVect, test, iCol, val;
    for fVect in Clique
    do
      test:=1;
      for iCol in [1..n]
      do
        val:=(eVect[iCol]-fVect[iCol]) mod 4;
        if val=2 then
          test:=0;
        fi;
        if val=3 then
          test:=0;
        fi;
      od;
      if test=1 then
        return true;
      fi;
    od;
    return false;
  end;


  for eVert in BuildSet(n, [0,1,2,3])
  do
    if IsCovered(eVert)=false then
      Add(ListVert, eVert);
    fi;
  od;
  return ListVert;
end;
