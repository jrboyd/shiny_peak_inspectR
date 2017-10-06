runx_bws = c("RUNX1" = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx1_pooled_FE.bw",
             "RUNX2" = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx2_pooled_FE.bw")

library(magrittr)
MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
  sub(pattern = "_FE.bw", replacement = "")
library(rtracklayer)
bw_gr = import.bw(MCF7bza_bws[1])
rng = range(bw_gr)
end(rng) = ceiling(end(rng)/100) * 100
win = slidingWindows(rng, width = 100, step = 100)
print(object.size(win), units = "GB")
print(object.size(bw_gr), units = "GB")
win = unlist(win)
sn = names(seqlengths(win))
names(sn) = sn
library(pbapply)
new_seqlengths = pbsapply(sn, function(x){
  max(end(subset(win, seqnames == x)))
})
seqlengths(win) = new_seqlengths
seqlengths(bw_gr) = seqlengths(win)
mid_gr = function(gr){
  start(gr) + floor((width(gr) - 1)/2)
}
mids = mid_gr(win)
start(win) = mids
end(win) = mids
olaps = findOverlaps(win, bw_gr)
win$FE = bw_gr[subjectHits(olaps)]$score
