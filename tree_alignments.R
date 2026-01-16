# setup -----------------------------------------------------------------------------------------------------------
library(tidyverse)
library(pegas)
library(here)
library(Biostrings)
library(msa)
library(fs)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)
library(pander)
library(paletteer)
library(patchwork)
library(ggfun)
library(seqinr)

# load('output/image.Rdata')

# stolen from shaunpwilkinson/insect
duplicated.DNAbin <- function(x, incomparables = FALSE, pointers = TRUE, ...){
  incomparables <- incomparables #TODO
  if(is.list(x)){
    if(pointers){
      hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...) # logical length x
      pntrs <- .point(hashes)
      attr(dupes, "pointers") <- pntrs
    }else{
      dupes <- duplicated(lapply(x, as.vector)) # as.vector removes attributes
    }
    names(dupes) <- names(x)
  }else if(is.matrix(x)){
    if(pointers){
      hashes <- apply(x, 1, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...)
      pntrs <- .point(hashes)
      attr(dupes, "pointers") <- pntrs
    }else{
      dupes <- duplicated(x, MARGIN = 1, ... = ...)
    }
    names(dupes) <- rownames(x)
  }else{
    dupes <- FALSE
    if(pointers) attr(dupes, "pointers") <- 1
  }
  return(dupes)
}

# stolen from shaunpwilkinson/insect
.point <- function(h){
  uh <- unique(h)
  pointers <- seq_along(uh)
  names(pointers) <- uh
  unname(pointers[h])
}

# stolen from shaunpwilkinson/insect
unique.DNAbin <- function(x, incomparables = FALSE, attrs = TRUE,
                          drop = FALSE, ...){
  incomparables <- incomparables
  if(is.list(x)){
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == length(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    dupes <- duplicated(lapply(x, as.vector), ... = ...)
    x <- x[!dupes]
    if(attrs) {
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!dupes]
      attributes(x) <- tmpattr
    }
  }else if(is.matrix(x)){
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == nrow(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    dupes <- duplicated(x, MARGIN = 1, ... = ...)
    x <- x[!dupes, , drop = drop]
    if(attrs){
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!dupes]
      attributes(x) <- tmpattr
    }
  }else{
    if(attr) x <- as.vector(x)
  }
  return(x)
}

# partition finder config template 
pft <- "
## ALIGNMENT FILE ##
alignment = {alignment_phy};

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = all;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
{data_block}

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;
"


# create all-species alignments -----------------------------------------------------------------------------------
rot <- function(x, n = 1) {
  n <- n %% length(x)             
  if (n == 0) return(x)
  c(tail(x, -n), head(x, n))
}

tips <- read_lines(here('data','valid_tips.txt'))

alignments %>%
  iwalk(\(seqs,marker) {
    seqs %>%
      iwalk(\(alignment,quality) {
        alignment <- alignment %>%
          unique()
        # tt <- tips[tips %in% labels(alignment)]
        # alignment <- alignment[tt,]
        samples <- all_samples %>%
          filter(id %in% rownames(alignment)) %>%
          group_by(species) %>%
          mutate( lbl = str_glue('{species}-{row_number()}')) %>%
          mutate( lbl = str_replace_all(lbl,'\\s+','_')) %>%
          arrange(lbl) %>%
          ungroup()
        alignment <- alignment[samples$id,]
        rownames(alignment) <- samples$lbl
        
        
        write_fasta(alignment,str_glue('output/alignments/{marker}_{quality}.fasta')) 
        write.nexus.data(alignment,str_glue('output/alignments/{marker}_{quality}.nex'),interleaved = FALSE) 
        write.dna(alignment,str_glue('output/alignments/{marker}_{quality}.phy'),format = "sequential",nbcol=-1,colsep="")
        
        frame <- 0:2 %>%
          map_int(\(frame) {
            f <- as.AAbin(apply(as.character(alignment),MARGIN = 1, FUN = \(seq) {
              seqinr::translate(seq,frame=frame,numcode=2)
            }))
            ifelse(sum(apply(f,MARGIN=1,FUN=\(row) sum(row == 42))) == 0,frame+1,NA)
          }) %>%
          keep(\(x) !is.na(x))
        
        cols <- ncol(alignment)
        start <- rot(1:3,frame-1)
        
        offset <- (cols - (frame-1)) %% 3
        end <- rot((cols-2):cols,-offset)
        
        data_block <- str_c(str_glue('codon_{1:3} = {start}-{end}\\3;'),collapse='\n')
        alignment_phy <- path_file(str_glue('output/alignments/{marker}_{quality}.phy'))
        
        pf_cfg <- str_glue(pft) 
        
        write_lines(pf_cfg,str_glue('output/alignments/{marker}_{quality}_partition_finder.cfg'))
      })
  })
# create goslinei-tapeinosoma concatenated alignment --------------------------------------------------------------
spp <- c('Plagiotremus goslinei','Plagiotremus tapeinosoma','Plagiotremus azaleus')
sss <- oll_samples %>%
  filter(species %in% spp)

sss_coi <- sss %>%
  filter(id %in% labels(alignments$coi$c50))
sss_cytb <- sss %>%
  filter(id %in% labels(alignments$cytb$c50))

unique_coi <- unique(alignments$coi$c50[sss_coi$id,])
coi_len <- ncol(unique_coi)
coi_missing <- str_c(rep('-',coi_len),collapse="")

unique_cytb <- unique(alignments$cytb$c50[sss_cytb$id,])
cytb_len <- ncol(unique_cytb)
cytb_missing <- str_c(rep('-',cytb_len),collapse="")

shared <- intersect(labels(unique_coi),labels(unique_cytb))
missing_from_cytb <- setdiff(labels(unique_coi),shared)
missing_from_coi <- setdiff(labels(unique_cytb),shared)

coi_missing_seqs <- coi_missing %>%
  rep(length(missing_from_coi)) %>%
  set_names(missing_from_coi) %>%
  DNAStringSet() %>%
  as.DNAbin() %>%
  as.matrix()

cytb_missing_seqs <- cytb_missing %>%
  rep(length(missing_from_cytb)) %>%
  set_names(missing_from_cytb) %>%
  DNAStringSet() %>%
  as.DNAbin() %>%
  as.matrix()

# write_fasta(cytb_missing_seqs,here('output','alignments','cytb_missing.fasta'))
# write_fasta(coi_missing_seqs,here('output','alignments','coi_missing.fasta'))

coi_full <- rbind(unique_coi,coi_missing_seqs)
coi_full <- coi_full[order(labels(coi_full)),]
nrow(coi_full)

cytb_full <- rbind(unique_cytb,cytb_missing_seqs)
cytb_full <- cytb_full[order(labels(cytb_full)),]
nrow(cytb_full)

all(labels(coi_full) == labels(cytb_full))

combined <- cbind(coi_full,cytb_full)
combined <- combined[order(labels(combined)),]

alignment_phy <- 'pgo-pta-paz-concatenated.phy'
alignment_fasta <- 'pgo-pta-paz-concatenated.fasta'
# write_fasta(combined,here('output','alignments',alignment_phy))
combined_samples <- sss %>%
  filter(id %in% labels(combined)) %>%
  drop_na(lon,lat) %>%
  st_as_sf(coords=c('lon','lat'),remove=FALSE,crs=4326) %>%
  st_join(regions,st_within) %>%
  mutate(meow_province = coalesce(meow_province,PROVINCE)) %>%
  mutate(meow_province = str_replace_all(meow_province,',','')) %>%
  dplyr::rename(name=id) %>%
  mutate(distribution = case_when(
    meow_province == 'Western Indian Ocean' ~ 'Western Indian Ocean',
meow_province == 'Red Sea and Gulf of Aden' ~ 'Red Sea',
meow_province == 'Hawaii' ~ 'Hawaii',
meow_province == 'Marquesas' ~ 'Marquesas',
REALM == 'Temperate Southern Africa' ~ 'Western Indian Ocean',
.default = REALM
  )) #%>%
# as_tibble()

combined_final <- combined[combined_samples$name,]
write.dna(combined_final,here('output','alignments',alignment_phy),format = "sequential",nbcol=-1,colsep="")
write_fasta(combined_final,here('output','alignments',alignment_fasta))
write_tsv(combined_samples %>% as_tibble() %>% select(name,distribution),here('output','alignments','pgo-pta-paz-concatenated_metadata.csv'))

coi_start <- 1
coi_end <- coi_len
cytb_start <- coi_len+1
cytb_end <- ncol(combined)

data_block <- str_trim(str_glue('
coi = {coi_start}-{coi_end}\\1;
cytb = {cytb_start}-{cytb_end}\\1;
'))

pfcfg <- str_glue(pft)
write_lines(pfcfg,here('output','alignments','pgo-pta-paz-concatenated_partition_finder.cfg'))


centroids <- combined_samples %>%
  group_by(distribution) %>%
  summarise(geometry = st_centroid(st_combine(geometry)))

distances <- centroids %>%
  st_distance()
rownames(distances) <- centroids$distribution
colnames(distances) <- centroids$distribution
write.csv(distances,here('output','alignments','pgo-pta-concatenated_distances.csv'),row.names = TRUE,quote = FALSE)

# alignment for goslinei/tapeinosoma beast analysis ---------------------------------------------------------------
spp <- c('Plagiotremus goslinei','Plagiotremus tapeinosoma','Plagiotremus azaleus')

alignment <- alignments %>%
  pluck('coi','c99')

samples <- all_samples %>%
  filter(species %in% spp) %>%
  filter(id %in% rownames(alignment)) %>%
  filter(!is.na(lat) & !is.na(lon))

alignment <- alignment[samples$id,] %>%
  unique()

samples <- samples %>%
  filter(id %in% rownames(alignment))

write_fasta(alignment,here('output','alignments','pgo-pta-beast.fasta'))
# alignment for ewaensis/rhinrhynchos beast analysis --------------------------------------------------------------
spp <- c('Plagiotremus ewaensis','Plagiotremus rhinorhynchos','Plagiotremus laudandus','Plagiotremus flavus')
# spp <- c('Plagiotremus ewaensis','Plagiotremus rhinorhynchos','Xiphasia matsubarai')
# spp <- c('Plagiotremus ewaensis','Plagiotremus rhinorhynchos','Omobranchus anolius')

alignment <- alignments %>%
  pluck('coi','c99')

samples <- all_samples %>%
  filter(species %in% spp) %>%
  filter(id %in% rownames(alignment)) %>%
  filter(!is.na(lat) & !is.na(lon))

alignment <- alignment[samples$id,] %>%
  unique()
samples <- samples %>%
  filter(id %in% rownames(alignment))

write_fasta(alignment,here('output','alignments','pew-prh-beast.fasta'))


# plagiotremus beast alignment ------------------------------------------------------------------------------------

# mixed haplotypes: XVII and III

spd <- all_samples %>%
  mutate(species_group = case_when(
    species == 'Plagiotremus tapeinosoma' & meow_province == 'Marquesas' ~ 'pta_marquesas',
    .default = species
  ))

unique_coi <- alignments$coi$c99 %>%
  unique()

spd <- spd %>%
  filter(id %in% rownames(unique_coi))

spda <- spd %>%
  group_by(species_group) %>%
  slice_sample(n = 4) %>%
  ungroup() %>%
  # filter(genus == 'Plagiotremus' | species == 'Xiphasia setifer')
  filter(genus == 'Plagiotremus' | species == 'Meiacanthus atrodorsalis')

beast_alignment <- unique_coi[spda$id,]
write_fasta(beast_alignment,here('output','alignments','plagiotremus-beast.fasta'))

# biogeography beast ----------------------------------------------------------------------------------------------

htt <- haplotypes$pgt$coi
alignment <- htt$alignment %>%
  unique()
samples <- htt$samples %>%
  filter(id %in% rownames(alignment)) %>%
  mutate(pvc = str_replace_all(province,'[[:space:]-]','')) %>%
  mutate(beast_id = str_glue("{id}_{pvc}"))

alignment <- alignment[samples$id,]
rownames(alignment) <- samples$beast_id

write_fasta(alignment,here('output','alignments','plagiotremus-beast_biogeography.fasta'))


# biogeography azaleus --------------------------------------------------------------------------------------------

# spp <- c('Plagiotremus goslinei','Plagiotremus tapeinosoma','Plagiotremus azaleus')
# 
# htt <- all_samples
# alignment <- alignments
# samples <- htt$samples %>%
#   filter(id %in% rownames(alignment)) %>%
#   mutate(pvc = str_replace_all(province,'[[:space:]-]','')) %>%
#   mutate(beast_id = str_glue("{id}_{pvc}"))
# 
# alignment <- alignment[samples$id,]
# rownames(alignment) <- samples$beast_id
# 
# write_fasta(alignment,here('output','alignments','plagiotremus-beast_biogeography.fasta'))

spp <- c('Plagiotremus goslinei','Plagiotremus tapeinosoma','Plagiotremus azaleus')

alignment <- alignments %>%
  pluck('coi','c99')

samples <- all_samples %>%
  filter(species %in% spp) %>%
  filter(id %in% rownames(alignment)) #%>%
  # filter(!is.na(lat) & !is.na(lon))

alignment <- alignment[samples$id,] %>%
  unique()

samples <- samples %>%
  filter(id %in% rownames(alignment)) %>%
  left_join(plagiotremus_data %>% select(id,province),by='id') %>%
  mutate(province = case_when(
    species == 'Plagiotremus azaleus' ~ 'Tropical East Pacific',
    .default = province
  )) %>%
  mutate(pvc = str_replace_all(province,'[[:space:]-]','')) %>%
  mutate(beast_id = str_glue("{id}_{pvc}"))

alignment <- alignment[samples$id,]
rownames(alignment) <- samples$beast_id

write_fasta(alignment,here('output','alignments','plagiotremus-beast_biogeography_azaleus.fasta'))
