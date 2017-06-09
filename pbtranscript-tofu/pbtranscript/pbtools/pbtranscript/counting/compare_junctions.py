import re

def overlaps(s1, s2):
    return max(0, min(s1.end, s2.end) - max(s1.start, s2.start))

def compare_junctions(r1, r2, group_info, fsm_maps, collapse_3_distance, collapse_5_distance, internal_fuzzy_max_dist=0):
    """
    r1, r2 should both be BioReaders.GMAPSAMRecord

    super
    exact
    subset
    partial
    nomatch

    <internal_fuzzy_max_dist> allows for very small amounts of diff between internal exons
    useful for chimeric & slightly bad mappings
    """

    # extract full-length group information
    g1, fl1, mfl1, fsm1 = 0, 0, -1, None
    for group in group_info[r1.seqid]:
        try:
            fsm1 = fsm_maps[group]
        except KeyError:
            pass

        aux1 = int( re.search( 'f.*p', group ).group(0)[1:-1] )
        #aux1 = int( re.search( 'f.*p', group.split( "|", 1 )[1] ).group(0)[1:-1] )
        if aux1 > mfl1:
            mfl1 = aux1
        fl1 += aux1 
        g1 += 1

    g2, fl2, mfl2, fsm2 = 0, 0, -1, None
    for group in group_info[r2.seqid]:
        try:
            fsm2 = fsm_maps[group]
        except KeyError:
            pass

        aux2 = int( re.search( 'f.*p', group ).group(0)[1:-1] )
        #aux2 = int( re.search( 'f.*p', group.split( "|", 1 )[1] ).group(0)[1:-1] )
        if aux2 > mfl2:
            mfl2 = aux2
        fl2 += aux2
        g2 += 1

    if fsm1 != None and fsm2 != None and fsm1 != fsm2:
        return "nomatch"

    #Set minimum distance to avoid collapses in 3' and 5'
    dist_l, dist_r = 0, 0
    if r1.strand == '+':
        dist_l, dist_r = collapse_5_distance, collapse_3_distance
    else:
        dist_l, dist_r = collapse_3_distance, collapse_5_distance

    # The same condition applied in compare exon matrix 
    if abs( r1.segments[0].start - r2.segments[0].start ) > dist_l:
        if r1.segments[0].start < r2.segments[0].start and fl2 > g2: 
            return "nomatch"
        if r1.segments[0].start > r2.segments[0].start and fl1 > g1:
            return "nomatch"

    if abs( r1.segments[-1].end - r2.segments[-1].end ) > dist_r:
        if r1.segments[-1].end < r2.segments[-1].end and fl1 > g1: 
            return "nomatch"
        if r1.segments[-1].end > r2.segments[-1].end and fl2 > g2:
            return "nomatch"

    found_overlap = False
    # super/partial --- i > 0, j = 0
    # exact/partial --- i = 0, j = 0
    # subset/partial --- i = 0, j > 0
    for i,x in enumerate(r1.segments):
        # find the first matching r2, which could be further downstream
        for j,y in enumerate(r2.segments):
            if i > 0 and j > 0: break
            if overlaps(x, y) > 0:
                found_overlap = True
                break
        if found_overlap: 
            break
    if not found_overlap: return "nomatch"
    # now we have r1[i] matched to r2[j]
    # if just one exon, then regardless of how much overlap there is, just call it exact
    if len(r1.segments) == 1:
        if len(r2.segments) == 1: return "exact"
        else:
            if r1.segments[0].end <= r2.segments[j].end:
                return "subset"
            else:
                return "partial"
    else:
        if len(r2.segments) == 1: return "super"
        else: # both r1 and r2 are multi-exon, check that all remaining junctions agree
            k = 0
            if ( i != 0 or j != 0 ):
                if abs(r1.segments[i+k].start-r2.segments[j+k].start)>internal_fuzzy_max_dist:
                    return "partial"
            while i+k+1 < len(r1.segments) and j+k+1 < len(r2.segments):
                if abs(r1.segments[i+k].end-r2.segments[j+k].end)>internal_fuzzy_max_dist or \
                   abs(r1.segments[i+k+1].start-r2.segments[j+k+1].start)>internal_fuzzy_max_dist:
                    return "partial"
                k += 1
            #print i, j, k
            if i+k+1 == len(r1.segments):
                if j+k+1 == len(r2.segments): 
                    if i == 0:
                        if j == 0: return "exact"
                        else: return "subset"    # j > 0
                    else: return "super"
                else: # r1 is at end, r2 not at end
                    if abs(r1.segments[i+k].end-r2.segments[j+k].end)>internal_fuzzy_max_dist:
                        return "partial"
                    if i == 0: return "subset"
                    else:  # i > 0
                        if abs(r1.segments[i+k-1].end-r2.segments[j+k-1].end)>internal_fuzzy_max_dist or \
                           abs(r1.segments[i+k].start-r2.segments[j+k].start)>internal_fuzzy_max_dist:
                            return "partial"
                        else: 
                            return "concordant"
            else: # r1 not at end, r2 must be at end
                if abs(r1.segments[i+k].end-r2.segments[j+k].end)>internal_fuzzy_max_dist:
                    return "partial"
                if j == 0: return "super"
                else:
                    if abs(r1.segments[i+k-1].end-r2.segments[j+k-1].end)>internal_fuzzy_max_dist or \
                        abs(r1.segments[i+k].start-r2.segments[j+k].start)>internal_fuzzy_max_dist:
                        return "partial"
                    else:
                        return "concordant"

