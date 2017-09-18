function inEntrezList = mapKeyVectToEntrezList(netobj,keyGeneVect)

    keyList = netobj.key(logical(keyGeneVect));
    inEntrezList = nanvalues(netobj.KeytoEntrezIDMap,keyList);

end