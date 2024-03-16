DihedralGroupPG:=function(r)
  return Group([PermList(List([1..r], x->NextIdx(r,x))), PermList(Reversed([1..r]))]);
end;

FuncComputeRelevantGroup:=function(r, eSet, Linner)
  local Dihedral, eStab, LinnerRed;
  Dihedral:=DihedralGroupPG(r);
  eStab:=Stabilizer(Dihedral, eSet, OnSets);
  LinnerRed:=Set(List(Linner, x->Set(x)));
  return Stabilizer(eStab, LinnerRed, OnSetsSets);
end;


ListBasicCases:=function(r)
  local ListSubset, i, Dihedral, O, ListCase, eO, eSet, eStab, eGrp, ListFinalCase, ListOpt, Linner, ListOptOrd, iNext, fSet, fStab, ListCaseInner, O2, eO2;
  ListSubset:=Combinations([1..r]);
  Dihedral:=DihedralGroupPG(r);
  O:=Orbits(Dihedral, ListSubset, OnSets);
  ListCase:=[];
  for eO in O
  do
    eSet:=eO[1];
    eStab:=Stabilizer(Dihedral, eSet, OnSets);
    ListOpt:=[];
    ListOptOrd:=[];
    for i in [1..r]
    do
      iNext:=NextIdx(r,i);
      if Position(eSet, i)<>fail and Position(eSet, iNext)<>fail then
        Add(ListOpt, [i, iNext]);
        Add(ListOptOrd, Set([i, iNext]));
      fi;
    od;
    ListCaseInner:=Combinations(Set(ListOptOrd));
    O2:=Orbits(eStab, ListCaseInner, OnSetsSets);
    for eO2 in O2
    do
      Linner:=List(eO2[1], x->ListOpt[Position(ListOptOrd, x)]);
      fStab:=Stabilizer(eStab, eO2[1], OnSetsSets);
      Add(ListCase, rec(eSet:=eSet, eBigGrp:=fStab, Linner:=Linner));
    od;
  od;
  return ListCase;
end;


RefinateListCase:=function(ListCase, n)
  local NewListCase, eCase, len;
  NewListCase:=[];
  for eCase in ListCase
  do
    len:=Length(eCase.eSet);
    if len>2 and len<n then
      Add(NewListCase, eCase);
    fi;
  od;
  return NewListCase;
end;



FuncGiveNextPair:=function(ListEdgeStatus, r, q, ePair)
  local NextPair, RML;
  if ePair[1]=PrevIdx(r, ePair[2]) then
    NextPair:=[ePair[2], NextIdx(r, ePair[2])];
  else
    NextPair:=[ePair[2], PrevIdx(r, ePair[2])];
  fi;
  RML:=First(ListEdgeStatus, x->x.iEdge=ePair[2]);
  if RML.eG="unset" then
    Print("You did something naughty, please be prudent\n");
    Print(NullMat(5));
  fi;
  return OnTuples(NextPair, RML.eG);
end;

IsDefinedNextPair:=function(ListEdgeStatus, ePair)
  return First(ListEdgeStatus, x->x.iEdge=ePair[2]).eG<>"unset";
end;




FuncGivePrevPair:=function(ListEdgeStatus, r, q, ePair)
  local PrevPair, RML;
  if ePair[2]=NextIdx(r, ePair[1]) then
    PrevPair:=[PrevIdx(r, ePair[1]), ePair[1]];
  else
    PrevPair:=[NextIdx(r, ePair[1]), ePair[1]];
  fi;
  RML:=First(ListEdgeStatus, x->x.iEdge=ePair[1]);
  if RML.eG="unset" then
    Print("You did something naughty, please be prudent\n");
    Print(NullMat(5));
  fi;
  return OnTuples(PrevPair, RML.eG);
end;

IsDefinedPrevPair:=function(ListEdgeStatus, ePair)
  return First(ListEdgeStatus, x->x.iEdge=ePair[1]).eG<>"unset";
end;

TheSaturation:=function(eGrp, ListEdgeStatus, TheSequence, r, q, IsInner, ListSetInner)
  local FuncMatch, NewSequence, ListPos, MaxPos, MinPos, MaxPair, MinPair, iMax, iMin, test, g, WeDidSomething, iPair, jPair, iPos, jPos, ePair1, ePair2, PrevEPair2, NextEPair1, iL;
  FuncMatch:=function(i)
    local U;
    for U in ListEdgeStatus
    do
      if U.iEdge=i then
        return true;
      fi;
    od;
    return false;
  end;
  NewSequence:=ShallowCopy(TheSequence);
  while(true)
  do
    ListPos:=List(NewSequence, x->x.Pos);
    MaxPos:=Maximum(ListPos);
    MinPos:=Minimum(ListPos);
    MaxPair:=First(NewSequence, x->x.Pos=MaxPos).ePair;
    MinPair:=First(NewSequence, x->x.Pos=MinPos).ePair;
    iMax:=MaxPair[2];
    iMin:=MinPair[1];
    # test that the inner status is respected
    if IsInner=true then
      if Position(ListSetInner, Set(MaxPair))=fail then
        return fail;
      fi;
    else
      if Position(ListSetInner, Set(MaxPair))<>fail then
        return fail;
      fi;
    fi;
    if IsInner=true then
      if Position(ListSetInner, Set(MinPair))=fail then
        return fail;
      fi;
    else
      if Position(ListSetInner, Set(MinPair))<>fail then
        return fail;
      fi;
    fi;
    # test that the edge status is correct
    if Length(NewSequence)>=2 then
      if FuncMatch(MaxPair[1])=false or FuncMatch(MinPair[2])=false then
        return fail;
      fi;
    fi;
    if IsInner=true then
      if FuncMatch(iMin)=false or FuncMatch(iMax)=false then
        return fail;
      fi;
    fi;
    # check pairwise consistency of the chain
    if Length(NewSequence)>=2 then
      for iL in [1..Length(NewSequence)-1]
      do
        ePair1:=First(NewSequence, x->x.Pos=MinPos+iL-1).ePair;
        ePair2:=First(NewSequence, x->x.Pos=MinPos+iL).ePair;
        if IsDefinedNextPair(ListEdgeStatus, ePair1) then
          NextEPair1:=FuncGiveNextPair(ListEdgeStatus, r, q, ePair1);
          g:=RepresentativeAction(eGrp, NextEPair1, ePair2, OnTuples);
          if g=fail then
            return fail;
          fi;
        fi;
        if IsDefinedPrevPair(ListEdgeStatus, ePair2) then
          PrevEPair2:=FuncGivePrevPair(ListEdgeStatus, r, q, ePair2);
          g:=RepresentativeAction(eGrp, PrevEPair2, ePair1, OnTuples);
          if g=fail then
            return fail;
          fi;
        fi;
      od;
    fi;
    # the length killing criterion
    test:=FuncMatch(iMax) or FuncMatch(iMin);
    if IsInner=false and Length(NewSequence)>q-1 then
      return fail;
    fi;
    if IsInner=false and Length(NewSequence)=q-1 and test=true then
      return fail;
    fi;
    if IsInner=false and Length(NewSequence)<=q-1 and test=false then
      return rec(Sequence:=NewSequence, IsFinished:=true);
    fi;
    if IsInner=true and Length(NewSequence)>=2*q+1 then
      for iPos in [MinPos..MaxPos-q]
      do
        jPos:=iPos+q;
        iPair:=First(NewSequence, x->x.Pos=iPos).ePair;
        jPair:=First(NewSequence, x->x.Pos=jPos).ePair;
        g:=RepresentativeAction(eGrp, iPair, jPair, OnTuples);
        if g=fail then
          return fail;
        fi;
      od;
      return rec(Sequence:=NewSequence, IsFinished:=true);
    fi;
    # the next one
    WeDidSomething:=false;
    if FuncMatch(iMax)=true and MaxPos<=q then
      if IsDefinedNextPair(ListEdgeStatus, MaxPair) then
        Add(NewSequence, rec(Pos:=MaxPos+1, ePair:=FuncGiveNextPair(ListEdgeStatus, r, q, MaxPair)));
        WeDidSomething:=true;
      fi;
    fi;
    if FuncMatch(iMin)=true and MinPos>=-q then
      if IsDefinedPrevPair(ListEdgeStatus, MinPair) then
        Add(NewSequence, rec(Pos:=MinPos-1, ePair:=FuncGivePrevPair(ListEdgeStatus, r, q, MinPair)));
        WeDidSomething:=true;
      fi;
    fi;
    if WeDidSomething=false then
      return rec(Sequence:=NewSequence, IsFinished:=false);
    fi;
  od;
end;

RepresentativesLeftCosets:=function(BigGrp, SmallGrp)
  return List(RightCosets(BigGrp, SmallGrp), x->Representative(x)^(-1));
end;

#
#
# We do the exhaustive enumeration of all
# cases of isohedral (r,q)-polycycles
EnumerationAllIsohedralCases:=function(eCase, r, q)
  local ListInner, ListSetInner, eGrp, Dihedral, ListEdgeStatus, iEdge, ListVertexSequence, iVert, IsInner, ListBefore, ListFinished, ListAfter, eCaseINC, eL, eCoset, eGchosen, NewListEdgeStatus, eStatus, jEdge, g, IsCorrect, NewListVertexSequence, IsFinished, eSequence, TheRES, TheNew, eElt;
  ListInner:=eCase.Linner;
  ListSetInner:=Set(List(ListInner, x->Set(x)));
  eGrp:=eCase.eGrp;
  Dihedral:=DihedralGroupPG(r);
  ListEdgeStatus:=List(eCase.eSet, x->rec(iEdge:=x, eG:="unset"));
  ListVertexSequence:=[];
  for iVert in [1..r]
  do
    if Position(ListInner, [iVert, NextIdx(r, iVert)])=fail then
      IsInner:=false;
    else
      IsInner:=true;
    fi;
    Add(ListVertexSequence, rec(IsInner:=IsInner, VLIST:=[rec(Pos:=0, ePair:=[iVert, NextIdx(r, iVert)])]));
  od;
  ListBefore:=[rec(ListEdgeStatus:=ListEdgeStatus, ListVertexSequence:=ListVertexSequence)];
  ListFinished:=[];
  while(true)
  do
    ListAfter:=[];
    for eCaseINC in ListBefore
    do
      eL:=First(eCaseINC.ListEdgeStatus, x->x.eG="unset");
      iEdge:=eL.iEdge;
      for eGchosen in RepresentativesLeftCosets(Dihedral, eGrp)
      do
        IsCorrect:=true;
        NewListEdgeStatus:=[];
        for eStatus in eCaseINC.ListEdgeStatus
        do
          jEdge:=eStatus.iEdge;
          g:=RepresentativeAction(eGrp, iEdge, jEdge, OnPoints);
          if g=fail then
            Add(NewListEdgeStatus, eStatus);
          else
            eElt:=g^(-1)*eGchosen;
            Add(NewListEdgeStatus, rec(iEdge:=jEdge, eG:=eElt));
            # if we have a plane of symmetry by the edge, then we should
            # have one on the other side as well.
            if Order(Stabilizer(eGrp, jEdge, OnPoints))<>Order(Stabilizer(eGrp, OnPoints(jEdge, eElt), OnPoints)) then
              IsCorrect:=false;
            fi;
          fi;
        od;
        NewListVertexSequence:=[];
        IsFinished:=true;
        for eSequence in eCaseINC.ListVertexSequence
        do
          TheRES:=TheSaturation(eGrp, NewListEdgeStatus, eSequence.VLIST, r, q, eSequence.IsInner, ListSetInner);
          if TheRES=fail then
            IsCorrect:=false;
          else
            Add(NewListVertexSequence, rec(IsInner:=eSequence.IsInner, VLIST:=TheRES.Sequence));
            if TheRES.IsFinished=false then
              IsFinished:=false;
            fi;
          fi;
        od;
        if IsCorrect=true then
          TheNew:=rec(ListEdgeStatus:=ShallowCopy(NewListEdgeStatus), ListVertexSequence:=NewListVertexSequence, eGrp:=eGrp);
          if IsFinished=true then
            Add(ListFinished, TheNew);
          else
            Add(ListAfter, TheNew);
          fi;
        fi;
      od;
    od;
    if Length(ListAfter)=0 then
      break;
    fi;
    ListBefore:=ShallowCopy(ListAfter);
  od;
  return ListFinished;
end;



FullEnumeration:=function(r, q)
  local LBAS, ListIsohedral, eCase, LIS;
  LBAS:=ListBasicCases(r);
  ListIsohedral:=[];
  for eCase in LBAS
  do
    if Length(eCase.eSet)>0 and Length(eCase.eSet)<r then
      Print("Working on eCase=", eCase, "\n");
      LIS:=EnumerationAllIsohedralCases(eCase, r, q);
      Append(ListIsohedral, LIS);
    fi;
  od;
  return ListIsohedral;
end;


IsSymmetryLiftable:=function(eCase, NewGroup)
  local ListEdgeStatus, ListiEdges, ListDone, i, iEdge, eGchosen, j, jEdge, RealChoice, ListAllChoices, g;
  ListEdgeStatus:=eCase.ListEdgeStatus;
  ListiEdges:=List(ListEdgeStatus, x->x.iEdge);
  ListDone:=ListWithIdenticalEntries(Length(ListiEdges), 0);
  for i in [1..Length(ListiEdges)]
  do
    iEdge:=ListiEdges[i];
    if ListDone[i]=0 then
      eGchosen:=ListEdgeStatus[i].eG;
      for j in [1..Length(ListiEdges)]
      do
        jEdge:=ListiEdges[j];
        g:=RepresentativeAction(NewGroup, iEdge, jEdge, OnPoints);
        if g<>fail then
          ListDone[j]:=1;
          RealChoice:=ListEdgeStatus[j].eG;
          ListAllChoices:=List(NewGroup, x->g^(-1)*eGchosen*x);
          if Position(ListAllChoices, RealChoice)=fail then
            return false;
          fi;
        fi;
      od;
    fi;
  od;
  return true;
end;




IsIsomorphicIsohedral:=function(eCase1, eCase2, BigGRP)
  local eGrp1, eGrp2, ConjGrp1, ListEdgeStatus1, ListEdgeStatus2, TestConjugacy, g;
  eGrp1:=eCase1.eGrp;
  eGrp2:=eCase2.eGrp;
  ListEdgeStatus1:=eCase1.ListEdgeStatus;
  ListEdgeStatus2:=eCase2.ListEdgeStatus;
  TestConjugacy:=function(g)
    local ListiEdges, ConjListEdgeStatus1, ListMultConj1, ListMult2, i;
    ListiEdges:=List(ListEdgeStatus2, x->x.iEdge);
    ConjListEdgeStatus1:=List(ListEdgeStatus1, x->rec(iEdge:=OnPoints(x.iEdge, g), eG:=g^(-1)*x.eG*g));
    ListMultConj1:=List(ListiEdges, x->First(ConjListEdgeStatus1, y->y.iEdge=x).eG);
    ListMult2:=List(ListiEdges, x->First(ListEdgeStatus2, y->y.iEdge=x).eG);
#    Print("ListMultConj1=", ListMultConj1, "\n");
#    Print("ListMult2=", ListMult2, "\n");
    for i in [1..Length(ListMult2)]
    do
      if RightCoset(eGrp2, ListMultConj1[i]^(-1))<>RightCoset(eGrp2, ListMult2[i]^(-1)) then
        return false;
      fi;
    od;
    return true;
  end;
  for g in BigGRP
  do
    ConjGrp1:=PersoGroupPerm(List(GeneratorsOfGroup(eGrp1), x->g^(-1)*x*g));
    if ConjGrp1=eGrp2 then
      if TestConjugacy(g)=true then
        return true;
      fi;
    fi;
  od;
  return false;
end;




AllGenerationFromAcase:=function(eCase, r, q)
  local eSet, eBigGrp, Linner, MaximalElementGroupFamily, ListCJ, ListCJrepr, ListSUB, eCJ, MaxSymmetryRepresentative, ListIsohedralMaximumSymmetry, FuncInsert, ListIso, MaxSym, fCase;
  MaximalElementGroupFamily:=function(TheFamily)
    local eSub, fSub, test;
    for eSub in TheFamily
    do
      test:=true;
      for fSub in TheFamily
      do
        if IsSubgroup(eSub, fSub)=false then
          test:=false;
        fi;
      od;
      if test=true then
        return eSub;
      fi;
    od;
    Print("We find nothing correct, something fishy\n");
    Print(NullMat(5));
  end;


  eSet:=eCase.eSet;
  eBigGrp:=eCase.eBigGrp;
  Linner:=eCase.Linner;
  ListCJ:=ConjugacyClassesSubgroups(eBigGrp);
  ListCJrepr:=List(ListCJ, x->Representative(x));
#  Print(NullMat(5));
  ListSUB:=[];
  for eCJ in ListCJ
  do
    Append(ListSUB, Set(eCJ));
  od;
  MaxSymmetryRepresentative:=function(eRep, GenSym)
    local ListSuperGroup, eSub, ListAdmisExtension;
    ListSuperGroup:=[];
    for eSub in ListSUB
    do
      if IsSubgroup(eSub, GenSym)=true then
        Add(ListSuperGroup, eSub);
      fi;
    od;
    ListAdmisExtension:=[];
    for eSub in ListSuperGroup
    do
      if IsSymmetryLiftable(eRep, eSub)=true then
        Add(ListAdmisExtension, eSub);
      fi;
    od;
    return MaximalElementGroupFamily(ListAdmisExtension);
  end;


  ListIsohedralMaximumSymmetry:=[];
  FuncInsert:=function(eCase, MaxSym)
    local TheInsert, eDone;
    TheInsert:=rec(ListEdgeStatus:=eCase.ListEdgeStatus, ListVertexSequence:=eCase.ListVertexSequence, eGrp:=MaxSym);
    for eDone in ListIsohedralMaximumSymmetry
    do
      if IsIsomorphicIsohedral(eDone, TheInsert, eBigGrp)=true then
        return;
      fi;
    od;
    Add(ListIsohedralMaximumSymmetry, TheInsert);
  end;

  for eCJ in ListCJrepr
  do
    ListIso:=EnumerationAllIsohedralCases(rec(eSet:=eSet, Linner:=Linner, eGrp:=eCJ), r, q);
    for fCase in ListIso
    do
      MaxSym:=MaxSymmetryRepresentative(fCase, eCJ);
      FuncInsert(fCase, MaxSym);
    od;
  od;
  return ListIsohedralMaximumSymmetry;
end;


DoingAndPresentingEnumeration:=function(r,q)
  local LBAS, LBAS2, LBAS3, ListFinished, eCase, AllGens, iCase, nbIsohedral;
  LBAS:=ListBasicCases(r);
  LBAS2:=Filtered(LBAS, x->Length(x.eSet)>0);
  LBAS3:=Filtered(LBAS2, x->Length(x.eSet)<r);
  Print("We have ", Length(LBAS3), " cases to consider\n");
  ListFinished:=[];
  for iCase in [1..Length(LBAS3)]
  do
    eCase:=LBAS3[iCase];
    Print("Treating iCase=", iCase, "\n");
    AllGens:=AllGenerationFromAcase(eCase, r, q);
    Add(ListFinished, rec(eSet:=eCase.eSet, Linner:=eCase.Linner, AllGens:=AllGens));
  od;
  nbIsohedral:=0;
  for iCase in [1..Length(LBAS3)]
  do
    eCase:=LBAS3[iCase];
    AllGens:=ListFinished[iCase].AllGens;
    if Length(AllGens)>0 then
      Print("i=", iCase, " eSet=", eCase.eSet, " Linner=", eCase.Linner, " |Sol|=", Length(AllGens), "\n");
      nbIsohedral:=nbIsohedral+Length(AllGens);
    fi;
  od;
  Print("Total of ", nbIsohedral, " isohedral (", r, ",", q, ")-polycycles\n");
  return ListFinished;
end;


MostDifficultTestingCase:=function(r,q)
  local LBAS, LBAS3, ListFinished, eCase, AllGens, iCase, nbIsohedral;
  LBAS:=ListBasicCases(r);
  LBAS3:=Filtered(LBAS, x->Length(x.eSet)=r);
  Print("We have ", Length(LBAS3), " cases to consider\n");
  ListFinished:=[];
  for iCase in [1..Length(LBAS3)]
  do
    eCase:=LBAS3[iCase];
    Print("Treating iCase=", iCase, "\n");
    AllGens:=AllGenerationFromAcase(eCase, r, q);
    Add(ListFinished, rec(eSet:=eCase.eSet, Linner:=eCase.Linner, AllGens:=AllGens));
  od;
  nbIsohedral:=0;
  for iCase in [1..Length(LBAS3)]
  do
    eCase:=LBAS3[iCase];
    AllGens:=ListFinished[iCase].AllGens;
    if Length(AllGens)>0 then
      Print("i=", iCase, " eSet=", eCase.eSet, " Linner=", eCase.Linner, " |Sol|=", Length(AllGens), "\n");
      nbIsohedral:=nbIsohedral+Length(AllGens);
    fi;
  od;
  if nbIsohedral<>1 then
    Print("There is a problem it seems\n");
  else
    Print("Things appear to be ok\n");
  fi;
end;



VertexOrbits:=function(r,q,eCase)
  local GRP, ListVertex, ListVertexOrd, GRA, eO, i, iElt, jElt, VertexSeq;
  GRP:=eCase.eGrp;
  ListVertex:=List([1..r], x->[x, NextIdx(r,x)]);
  ListVertexOrd:=List(ListVertex, x->Set(x));
  GRA:=NullGraph(Group(()), r);
  for eO in Orbits(GRP, ListVertexOrd, OnSets)
  do
    for i in [1..Length(eO)-1]
    do
      iElt:=Position(ListVertexOrd, eO[i]);
      jElt:=Position(ListVertexOrd, eO[i+1]);
      if iElt<>jElt then
        AddEdgeOrbit(GRA, [iElt, jElt]);
        AddEdgeOrbit(GRA, [jElt, iElt]);
      fi;
    od;
  od;
  for VertexSeq in eCase.ListVertexSequence
  do
    for i in [1..Length(VertexSeq.VLIST)-1]
    do
      iElt:=Position(ListVertexOrd, Set(VertexSeq.VLIST[i].ePair));
      jElt:=Position(ListVertexOrd, Set(VertexSeq.VLIST[i+1].ePair));
      if iElt<>jElt then
        AddEdgeOrbit(GRA, [iElt, jElt]);
        AddEdgeOrbit(GRA, [jElt, iElt]);
      fi;
    od;
  od;
  return List(ConnectedComponents(GRA), x->ListVertex{x});
end;



IsogonalCases:=function(r,q)
  local LBAS, i, eSet, Linner, j, eBigGrp, ListFinished, iCase, eCase, AllGens, nbIsogonal;
  LBAS:=[];
  for i in [2..q]
  do
    eSet:=[1..i];
    Linner:=[];
    for j in [1..i-1]
    do
      Add(Linner, [j, j+1]);
    od;
    eBigGrp:=FuncComputeRelevantGroup(q+1, eSet, Linner);
    Add(LBAS, rec(eSet:=eSet, Linner:=Linner, eBigGrp:=eBigGrp));
  od;
  ListFinished:=[];
  for iCase in [1..Length(LBAS)]
  do
    eCase:=LBAS[iCase];
    Print("Treating iCase=", iCase, "\n");
    AllGens:=AllGenerationFromAcase(eCase, q+1, r);
    Add(ListFinished, rec(eSet:=eCase.eSet, Linner:=eCase.Linner, AllGens:=AllGens));
  od;
  nbIsogonal:=0;
  for iCase in [1..Length(LBAS)]
  do
    eCase:=LBAS[iCase];
    AllGens:=ListFinished[iCase].AllGens;
    if Length(AllGens)>0 then
      Print("i=", iCase, " eSet=", eCase.eSet, " Linner=", eCase.Linner, " |Sol|=", Length(AllGens), "\n");
      nbIsogonal:=nbIsogonal+Length(AllGens);
    fi;
  od;
  Print("Total of ", nbIsogonal, " isogonal (", r, ",", q, ")-polycycles\n");
  return ListFinished;
end;


SearchIsohedralIsogonal:=function(r,q)
  local LBAS, LBAS2, LBAS3, ListFinished, eCase, AllGens, iCase, nbIsoHedralGonal, ListIsogonal;
  LBAS:=ListBasicCases(r);
  LBAS2:=Filtered(LBAS, x->Length(x.eSet)>0 and Length(x.eSet)<r and Length(x.Linner)=0);
  Print("We have ", Length(LBAS2), " cases to consider\n");
  ListFinished:=[];
  for iCase in [1..Length(LBAS2)]
  do
    eCase:=LBAS2[iCase];
    Print("Treating iCase=", iCase, "\n");
    AllGens:=AllGenerationFromAcase(eCase, r, q);
    ListIsogonal:=Filtered(AllGens, x->Length(VertexOrbits(r,q, x))=1);
    Add(ListFinished, rec(eSet:=eCase.eSet, Linner:=eCase.Linner, AllGens:=ListIsogonal));
  od;
  nbIsoHedralGonal:=0;
  for iCase in [1..Length(LBAS2)]
  do
    eCase:=LBAS2[iCase];
    AllGens:=ListFinished[iCase].AllGens;
    if Length(AllGens)>0 then
      Print("i=", iCase, " eSet=", eCase.eSet, " Linner=", eCase.Linner, " |Sol|=", Length(AllGens), "\n");
      nbIsoHedralGonal:=nbIsoHedralGonal+Length(AllGens);
    fi;
  od;
  Print("Total of ", nbIsoHedralGonal, " isohedral and isogonal (", r, ",", q, ")-polycycles\n");
  return ListFinished;
end;


ThreeNumbers:=function(r,q)
  local LS1, LS2, LS3, nb1, nb2, nb3, FuncNumber;
  FuncNumber:=function(LS)
    local TheNb, eCase;
    TheNb:=0;
    for eCase in LS
    do
      TheNb:=TheNb+Length(eCase.AllGens);
    od;
    return TheNb;
  end;
  LS1:=DoingAndPresentingEnumeration(r,q);
  LS2:=IsogonalCases(r,q);
  LS3:=SearchIsohedralIsogonal(r,q);
  nb1:=FuncNumber(LS1);
  nb2:=FuncNumber(LS2)-1; # because we exclude r-gon
  nb3:=FuncNumber(LS3);
  Print("[", nb1, ",", nb2, ",", nb3, "]");
  Print("  isohedral, isogonal, isohedral and isogonal (", r, ",", q, ")-polycycles\n");
end;



#
# We have some detected isohedral isogonal (r,q)-polycycle,
# which have vertices of different vertex degree.
# We need to hunt for the bug.
BUGSEARCH_SearchIsohedralIsogonal:=function(r,q)
  local LBAS, LBAS2, LBAS3, ListFinished, eCase, AllGens, iCase, nbIsoHedralGonal, ListIsogonal, ListBaddies;
  LBAS:=ListBasicCases(r);
  LBAS2:=Filtered(LBAS, x->Length(x.eSet)>0 and Length(x.eSet)<r and Length(x.Linner)=0);
  Print("We have ", Length(LBAS2), " cases to consider\n");
  ListFinished:=[];
  for iCase in [1..Length(LBAS2)]
  do
    eCase:=LBAS2[iCase];
    Print("Treating iCase=", iCase, "\n");
    AllGens:=AllGenerationFromAcase(eCase, r, q);
    ListIsogonal:=Filtered(AllGens, x->Length(VertexOrbits(r,q, x))=1);
    ListBaddies:=Filtered(ListIsogonal, x->Length(Set(List(x.ListVertexSequence, y->Length(y.VLIST))))<>1);
    if Length(ListBaddies)>0 then
      Print("We got one!!!!\n");
      Add(ListFinished, rec(eSet:=eCase.eSet, Linner:=eCase.Linner, AllGens:=ListBaddies));
    fi;
  od;
  return ListFinished;
end;
