function  [endo, epi, aod1, aod2] = adjust_position(endo, epi, aod1, aod2)

    if(length(endo)==1)
        endop = endo.vertices;
        epip = epi.vertices;
    else
        endop = endo;
        epip = epi;
    end
    added = 0;

    diff = max(epip(:,3)) - max(endop(:,3)) + added;
    endop(:,3) = endop(:,3) + repmat(diff,length(endop),1);
    aod1(:,3) = aod1(:,3) + repmat(diff,length(aod1),1);
    
    endo = endop;
    

end