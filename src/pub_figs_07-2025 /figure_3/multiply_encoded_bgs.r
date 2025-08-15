require(foreach)
require(tidyverse)
require(seqinr)
require(adagio)
require(cowplot)
require(viridis)
require(RColorBrewer)
require(extrafont)
#install.packages('lpSolve', repos='http://cran.rstudio.com')
#install.packages('adagio')

#setwd("~/Dropbox/active/hiv_bnabs_parallelism/code/")

#This file is the multiple encoding analysis run for figure 3 panel a

indir = '../../../data/clyde_westfall_2024_final/'
outdir = '../../../results/pub_figs_03-2025/figure_3/multi_encode/'

#this is from Zanini 2017
mutations <- bind_rows(
        tibble(nucmut = "g>c", rate = 10^-7, type = "tv"), 
        tibble(nucmut = "a>t", rate = 7 * 10^-7, type = "tv"), 
        tibble(nucmut = "c>g", rate = 5 * 10^-7, type = "tv"), 
        tibble(nucmut = "a>c", rate = 9 * 10^-7, type = "tv"), 
        tibble(nucmut = "g>t", rate = 2 * 10^-6, type = "tv"), 
        tibble(nucmut = "t>a", rate = 3 * 10^-6, type = "tv"), 
        tibble(nucmut = "t>g", rate = 3 * 10^-6, type = "tv"), 
        tibble(nucmut = "c>a", rate = 5 * 10^-6, type = "tv"), 
        tibble(nucmut = "a>g", rate = 6 * 10^-6, type = "ts"), 
        tibble(nucmut = "t>c", rate = 10^-5, type = "ts"), 
        tibble(nucmut = "c>t", rate = 1.2 * 10^-5, type = "ts"), 
        tibble(nucmut = "g>a", rate = 1.6 * 10^-5, type = "ts"))

#this will plot some sanity checks in the loop if TRUE
debug.bool <- TRUE

runParams <- list()
runParams[['3BNC']] <- list(indir = paste0(indir, '3BNC117/'), 
                            regexp = "^2[CDE][0-9]$", 
                            save = "3BNC117")
runParams[['10-1074']] <- list(indir = paste0(indir, '10-1074/'), 
                            regexp = "^1H[BCD]?[0-12]+[K]?", 
                            save = "10-1074")



#for both 10-1074 and 3BNC117
foreach(bnab = names(runParams))%do%{
    print('started foreach loop again :/')  
  
    indir <- runParams[[bnab]]$indir
    regexp <- runParams[[bnab]]$regexp
    saveval <- runParams[[bnab]]$save
    
    pts <- list.files(indir)
    pts <- pts[grep(regexp, pts)]

    encodings <- foreach(pt = pts, .combine = "rbind")%do%{

        infa <- list.files(paste0(indir, pt))
        #Get the untranslated FASTA file
        fa.fi <- infa[grep("translated", infa, invert = TRUE)]

        fa <- read.fasta(paste0(indir, pt, "/", fa.fi))

        #Translate the fasta files in frame and put the results in a tibble
        #Each row is a sequence and each column is an AA position
        AApos <- as_tibble(t(as_tibble(lapply(fa, translate)))) %>% 
            #Add back in the sequence identifiers 
            mutate(seq = names(fa)) 

        #Similar as above, but for nucleotides as opposed to amino acids
        nucpos <- as_tibble(t(as_tibble(fa))) %>% 
            mutate(seq = names(fa)) 

        #Convert wide to long format for AA data and make position numeric
        aa_inf <- AApos %>% gather(ind, AA, -seq) %>% 
            mutate(AA_pos = as.numeric(gsub("V", "", ind))) %>%
            select(-ind)
        #Note - this is the ARRAY index, NOT hxb2 indexing at this point
        
        #This is a stupid way to get the nucleotide length
        maxind <- max((nucpos %>% gather(ind, nuc, -seq) %>% 
                       mutate(ind = as.numeric(gsub("V", "", ind))))$ind)

        #this is a smarter way to get the nucleotide length but I haven't tested that 
        #they're always the same
        maxind <- ncol(nucpos) - 1

        #I want the largest multiple of 3 that is less than the nucleotide length
        newmaxind <- floor(maxind/3)*3

        if(debug.bool == TRUE){

            #visually inspect that sequences are aligned
            #Note - this is the ARRAY index, NOT hxb2 indexing at this point
            aa_inf %>% filter(between(AA_pos, 332, 334)) %>%  
                ggplot() + geom_tile(aes(y = seq, x = AA_pos, fill = AA))
        }
        print("made it past data read in")

#Ok, in 3BNC, I'm running into a problem where we have non-multiples of three

        nuc_inf <- nucpos %>% 
            #convert nucpos (nucleotide position) to long
            gather(ind, nuc, -seq) %>% 
            mutate(ind = as.numeric(gsub("V", "", ind))) %>% 
            #and then filter anything that is greater than the final codon
            filter(ind <= newmaxind) %>% 
            #Rearrange the results to cluster by sequence and sort 
            group_by(seq) %>% arrange(seq, ind) %>% 
            #add information about which AA position each nuc belongs to (AA_pos)
            #and which position within the codon (pos))
            mutate(AA_pos = sort(rep(1:(n()/3), 3)), pos = rep(1:3, n()/3)) %>% 
            #drop the former nucleotide array index
            select(-ind)

                  
        #merge each row to contain info about each nuc and each AA
        full_inf <- full_join(aa_inf, nuc_inf) %>% 
            #grab the day/week number and retain only the first character
            #the reason that it's ok to only keep one char is that we only
            #need to check if it's day 0 or not
            mutate(time = substr(gsub("^.*_[d|D|w|W]", "", seq), 0, 1)) %>% 
            #filter out the hxb2 sequence
            filter(time != "H") %>% 
            #and convert time to numeric
            mutate(time = as.numeric(time))
    
#Lots of ways this could go wrong... will need to check these
        
        #Here, we're looking for the AA pattern right at the glycosolation site
        #because we're going to want to compare 332, 334, 325 across participants
        #we need a common indexing scheme

        #The occurs function is looking for the location of a sequence in 
        #a longer array. However, it doesn't work on characters, so we convert both
        #HXB2 and the alphabet to their numeric equivalents
        cnis <- occurs(match(c("C", "N", "I", "S"), LETTERS),
                       match(translate(fa[[1]]), LETTERS))
        #cnis is the AA array position of the C

        #make a tibble dictionary of the escape sites and their array positions
        pos_inf <- tibble(AA_pos = c(cnis[1] - 6, cnis[1]+1, cnis[1] + 3), 
                          ID = c(325, 332, 334))
        
        #convert full_inf back to wide format so that each nuc position
        #is in the same row of a single AA position
        semi <- full_inf %>% spread( pos, nuc) %>% 
            mutate(triplet = paste0(`1`, `2`, `3`))  %>%
            arrange(seq, AA_pos)
        #note, this is named semi because it's intermediate between wide and long


#quick vis - does this look right?
        if(debug.bool == TRUE){
            semi %>% filter(between(AA_pos, 325, 334)) %>%  
                ggplot() + geom_tile(aes(y = seq, x = AA_pos, fill = AA))
        }

    
    #I don't want anything intermediate frequency? Or maybe I do, I don't know
        
        #Here, we are trying to determine what the ancestral versions of 
        #each AA was at time 0 in the trial
        mono <- semi %>% filter(time == 0) %>% filter(AA != "X") %>% 
            #count the number of occurrences of each codon and each site
            group_by(AA_pos, AA, triplet) %>% dplyr::count() %>% 
            #group at the site level
            group_by(AA_pos) %>% 
            #count the largest category (catn), 
            #the total number of sequences with an identity at that site (N)
            #the highest frequency category (topf)
            #and the frequency for each triplet encoding (truef)
            mutate(catn = max(n), N = sum(n), topf = catn/N, truef = n/N)  %>% 
            #retain only a few of those categories (do we need topf??)
            select(AA_pos, AA, triplet, truef) %>% 
            rename(ancAA = AA, ancTriplet = triplet)

           
        #EVR/AFF stopped code reviewing here on Jan 29, 2025

        if(debug.bool == TRUE){
            #plot AA frequencies by position
            mono %>% filter(AA_pos < 20) %>% ggplot() + geom_point(aes(x = AA_pos, y = truef))
            #all of the encodings at a given position sum to 1
            mono %>% group_by(AA_pos) %>% summarize(truef = sum(truef)) %>% 
                ggplot() + geom_point(aes(x = AA_pos, y = truef))

        }

    #Now let's look at derived alleles

        #this is maybe deletable - it's not saved, and doesn't seem 
        #to be exactly what I want
        left_join(semi %>% filter(time > 0), mono, 
                  relationship = "many-to-many") %>% 
            filter(!is.na(ancAA), AA != "X") %>% filter(AA != ancAA)

        filter_freq <- 0

        #merge together each position's info for derived and anc AAs
        #separate by sequence. truef is the day 0 frequency across seqs.
        derived_intermediate <- left_join(semi %>% filter(time > 0), 
                             mono %>% filter(truef >= filter_freq),
                             relationship = "many-to-many") %>% 
            filter(!is.na(ancAA), AA != "X")

        #We're splitting the ancestral AAs into nucs, and then comparing
        #by codon position, whether the anc nuc and current nuc differ
        derived <- derived_intermediate   %>%
            separate(ancTriplet, into = c('anc0', 'anc1', 'anc2', 'anc3'), 
                     sep = "", remove = FALSE) %>% 
            mutate(pos1_muts = ifelse(`1` != anc1, paste0(anc1, ">", `1`), ""), 
                   pos2_muts = ifelse(`2` != anc2, paste0(anc2, ">", `2`), ""), 
                   pos3_muts = ifelse(`3` != anc3, paste0(anc3, ">", `3`), ""),
                   AAmut = paste0(ancAA, ">", AA), 
                   nucmut = paste0(pos1_muts, pos2_muts, pos3_muts)) %>% 
            select(seq, AA_pos, time, ancAA, AA, AAmut, ancTriplet, triplet, 
                   nucmut, truef) 


        derived_new <- derived %>% #filter(AA_pos == 339) %>% 
            #Retain only the most closely related ancestor to the current seq
            #note, we think it's not necessary to group by seq here - same ans
            group_by(AA_pos, triplet) %>% 
            filter(nchar(nucmut) == min(nchar(nucmut))) %>% 
            #in the case of equidistant ancestors, retain whichever was most
            #common at day 0
            filter(truef == max(truef)) %>% 
            #remove sites that didn't mutate
            filter(triplet != ancTriplet) %>% ungroup() %>% 
            #remove sites that are more than 1 mutation away from an ancestor
            filter(nchar(nucmut) == 3) %>% 
            #remove synonymous sites 
            filter(ancAA != AA) 

        #EVR/AFF stopped code reviewing here on Feb 19, 2025

       
    #I want to know how many derived alleles are only one mutation
#away from multiple ancestral haplotypes

#Ok, let's say that in the case of equidistance, let's 
#choose the more frequent haplotype

        if(debug.bool == TRUE){

    #These are overall abundances
            g1 <- left_join(derived %>% filter(nchar(nucmut) == 3), mutations) %>% 
                ggplot() + geom_bar(aes(x = nucmut, fill = rate))
            

    #this is only counting each site once per mutation
            g2 <- left_join(derived %>% filter(nchar(nucmut) == 3), mutations) %>%
                group_by(AA_pos, nucmut, rate) %>% dplyr::count() %>% ungroup() %>% 
                ggplot() + geom_bar(aes(x = nucmut, fill = rate))

    #nice
            plot_grid(g1, g2)

        }


    #How many of these are multiply encoded? 
        #note, removed redundant filtering condition

        #Add the mutation rates and the equivalent HXB2 coordinates as "ID"
        #only for positions 325, 332, 334 - left_join will auto-fill NAs for 
        #other sites
        withmu <- left_join(left_join(derived_new, 
                                      mutations), pos_inf)

        if(debug.bool == TRUE){

            withmu %>% group_by(AA_pos, AAmut, nucmut, triplet, rate, ID) %>% 
                dplyr::count() %>% group_by(AA_pos, AAmut, ID) %>% mutate(N = sum(n)) %>% 
                    filter(N > 1) %>% 
                    ggplot() + 
                    geom_bar(aes(x = paste0(AA_pos,  "_", AAmut, "_", ID), y = n, fill = rate), 
                             stat = "identity") + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

        }

        #Label the participant from which these sequences are from and return
        withmu %>% mutate(pt = pt)

    }

    #save the table
    print('Saving the table')
    print(paste0(outdir, saveval,"_encodings.csv"))
    write.table(encodings %>% mutate(type = saveval), 
            file = paste0(outdir, saveval,"_encodings.csv"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
    print('Saved the table')

}
print('Made it out of loop')




#there are really not very many multiply encodable sites. 
#here, we are exhaustively searching all codon to codon transitions for those that are
#separated by a single nucleotide change, so we can see what AA transitions are 
#multiply encoded (i.e., AAT/AAC >AAG/AAA is K>N??*double check this)

nucs <- c("a", "t", "c", "g")
allnucs <- tibble(expand.grid(nucs, nucs, nucs, nucs, nucs, nucs)) 
names(allnucs) <- c('a1', 'a2', 'a3', 'd1', 'd2', 'd3')


allmutopts <- 
#Filter down to triplets separated by exactly one nucleotide mutation
allnucs %>% mutate(difs = (a1 != d1) + (a2 != d2) + (a3 != d3)) %>% 
    filter(difs == 1) %>% 
#Record the exact the nucleotide change in a column corresponding to its codon position
    mutate(ancTriplet = paste0(a1, a2, a3), triplet = paste0(d1, d2, d3)) %>% 
    mutate(pos1_muts = ifelse(a1 != d1, paste0(a1, ">", d1), ""), 
           pos2_muts = ifelse(a2 != d2, paste0(a2, ">", d2), ""), 
           pos3_muts = ifelse(a3 != d3, paste0(a3, ">", d3), ""),
           #and then also in a column aggregating across all codon positions
           nucmut = paste0(pos1_muts, pos2_muts, pos3_muts)) %>% 
#and also record what codon position was mutated
    mutate(mutpos = ifelse(a1 != d1, 1, 
                    ifelse(a2 != d2, 2, 
                    ifelse(a3 != d3, 3, NA)))) %>% 
#separately for each row, translate the nucleotides to AAs
    rowwise() %>% 
    mutate(ancAA = translate(c(as.character(a1), as.character(a2), as.character(a3))),
           AA = translate(c(as.character(d1), as.character(d2), as.character(d3)))) %>%
    select(ancTriplet, triplet, nucmut, mutpos, ancAA, AA) %>% 
#and then only retain things that are non-syn changes
    filter(ancAA != AA) %>% mutate(AAmut = paste0(ancAA, ">", AA))
   
multis <- allmutopts %>% group_by(ancTriplet, AAmut) %>% 
    summarize(possible_enc = n()) 
#multis records, hypothetically, from a single ancestral codon, 
#what are the AA mutational changes that are possible, and via 
#how many nucleotide changes



#10-1074 analysis
    
#Let's read in that encoding data from earlier

print("reading the table")
encodings <- tibble(read.table(paste0(outdir, "10-1074_encodings.csv"), header = TRUE))

encoding_types <- encodings %>% 
    group_by(ancTriplet, AA_pos, pt, triplet, ID, rate, AAmut, nucmut) %>% 
    dplyr::count() %>% ungroup()

#encoding_types summarizes across sequences to get counts of how many times each 
#codon to codon transition occurs at a given site in a given participant across all post-treatment
#times

#Ok, the core information that we want for our 10-1074 plot is only at a few ancestral AAs 
saved_results <- left_join(encoding_types, multis) %>% filter(possible_enc > 1)
write.table(saved_results, 
            file = "multiple_encodings.csv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)



#First, we retain only ancestral codons that CAN be multiply encoded (as calculated in multis)
fig_101074_maintext <- left_join(encoding_types, multis) %>% 
    filter(possible_enc > 1) %>% 
    group_by(ancTriplet, AAmut) %>% 
    #then we manually filter down to those codons that occur at 325, 332, 334 (i.e., escape sites)
    filter(is.element(ancTriplet, c("gac", "gat", "aac", "aat", "agt"))) %>% 
    #Let's make a custom facet label so that we can order how we like
    mutate(facet = factor(paste0(ancTriplet,  "_", AAmut, "_", ID), 
                          levels = c("gat_D>E_325", "gat_D>E_NA", 
                                     "gac_D>E_325", "gac_D>E_NA", 
                                     #huh, looks like some of the N>Ks are at 325... 
                                     #that seems wrong- nope, it's not
                                     "aat_N>K_325", "aat_N>K_NA",
                                     "aac_N>K_332", "aac_N>K_NA", 
                                     "agt_S>R_334", "agt_S>R_NA"))) %>%
    ggplot() + 
    geom_bar(aes(x = paste0(AA_pos,  "_", pt, "_", ID), y = n, fill = triplet), 
             col = "white", 
             stat = "identity") + 
        facet_grid(~facet, 
                   scales = "free_x", space = "free") + 
    theme_classic() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    labs(y = "Count (n)", fill = element_blank())

theme_AC <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_line(color="lightgrey", size=0.2),
      panel.grid.major.y = element_line(color="lightgrey", size=0.2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=18, family=font),
      axis.title = element_text(size=18, family=font),
      #legend.title=element_text(size=16, family=font, face="bold"),
      #legend.position="right",
      plot.title = element_text(size=22, family=font)
    )
}

#extrafont::font_import()

extrafont::loadfonts()

#I'm also generating a version that has as much as possible stripped away from the plot
#so that it's easier to make custom legends in a graphics program
fig_101074_maintext_minimal <- fig_101074_maintext  + 
    #this cuts out a little extra white space at the bottom of the plot
    coord_cartesian(ylim = c(1, 40)) + 
    #Hmm, I'm not sure I actually like Abby's theme here, but you can mess around with what you like
    #theme_AC() + 
    theme(strip.text = element_blank()) + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    theme(
      axis.text = element_text(size=6),
      axis.title = element_text(size=6),
      legend.text = element_text(size=6)
    )+ 
  
    scale_fill_brewer(palette = 'Set2') +
    theme(legend.key.size = unit(0.2, "cm"))+ 
    theme(panel.spacing = unit(0.1, "lines"))

    
#Honestly, for whatever reason, trying to explicitly specify Arial is not playing
#well with ggsave but specifying nothing produces Arial...
    ## theme(
    ##   axis.text = element_text(size=6, family="Arial"),
    ##   axis.title = element_text(size=6, family="Arial"),
    ##   plot.title = element_text(size=22, family="Arial")
    ## )


#full version with unpretty labels
ggsave(paste0(outdir, "DRM_encodings_10-1074.pdf"), fig_101074_maintext, 
       dpi = 300, width = 10, height = 2)

#version with no labels
#Labels were created manually
ggsave(paste0(outdir, "DRM_encodings_10-1074_minimal.pdf"), 
       fig_101074_maintext_minimal, dpi = 300, width = 5, height = 1.2)
embedFonts(file="../results/figs/DRM_encodings_10-1074_minimal.pdf")



#Now let's look at 3BNC - 
#We don't make a main text figure describing this, but we do pull out particular positions 
#that seem to be multiply encoded and put them in main text Figure 4

encodings <- tibble(read.table("../dat/3BNC117_encodings.csv", header = TRUE))

encoding_types <- encodings %>% 
    group_by(ancTriplet, AA_pos, pt, triplet, ID, rate, AAmut, nucmut) %>% 
    dplyr::count() %>% ungroup()

#Pull out all multiply encoded alleles
long <- encoding_types %>% group_by(AAmut, AA_pos, pt) %>% 
    mutate(numTypes = n()) %>% filter(n() > 1) %>% 
    arrange(pt, AA_pos, AAmut) 

#Convert from long to wide
wide <- long %>% select(ancTriplet, AA_pos, pt, AAmut, triplet, n) %>% 
    #rename derived triplet for clarity
    rename(derTriplet = triplet) %>% 
    group_by(AA_pos, pt) %>% arrange(desc(n)) %>% 
    mutate(derived = paste0("g", 1:n())) %>% 
    pivot_wider(names_from = derived, values_from = c(ancTriplet, derTriplet, n))

#We are also paying attention to cases in which we start with a CAA/CAG ancestral state
#and we go to a AAA/AAG derived state. 

#Let's write both wide and long formats of the data
write.table(long, file = paste0(outdir, "multiply_encoded_3BNC117_long.tsv"), 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(wide, file = paste0(outdir, "multiply_encoded_3BNC117_wide.tsv"), 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
