## make TSDA quercus tree
## ahipp@mortonarb.org, 2017-11-14

library(ape)
library(phytools)
library(ggtree)

# intro panel: source of the tree, subject to change,
# scope of research, what the tree does and doesn't tell you,
# importance of phylogenetic history

## Targeting lay-affinitive and industry audiences
## aim for 70-85 spp
## -- to add: insignis,
## aim for 4-5 stories
## ring around the edge: geographic regions (3-4 bands)
## we can do 4 bands of 1" tile
## think of the center sign and the nodes as a chance to interpret the geographic regions
## 3 curved benches, one for each clade: Quercus, Lobatae, OW oak clade
## vertical panels can have a leaf, an acorn, some other graphic in addition to the name:
##   -- will these be trees that you can point to or go see?
##   -- aim for 7-10 vertical panels

## AH TO DO:
##  1. cut down Ilex and Cyclobalanopsis a bit
##  2. Add in Mexican red and white oaks
##  3. rejigger signage

# Questions:
#  1. what information can we put on ind'l spp panels?
#    -- leaf outlines will be confused
#    -- shoot for 5 nodes plus a central panel
#    -- The central panel might have branches beneath it,
#       as it overlaps the base of the tree
#  2. Age ranges, different stories for different audiences?
#     -- still in the works
#  3. Trees around the edge, one per cUIClade? (at least 5 trees)
#    -- we have room to plant 4 more trees
#  4. paths wide enough for kids to walk?
#    -- these are taped on brick
#  ADDITIONAL NOTES FROM CONVERSATION WITH DAN GEORGE, 2017-11-14
#    -- intro and exit signs? people will walk both directions
#    -- sign on the north side comes out of "Trees in a changing climate"
#       ... tie into this, providing a deep-time perspective.
#       ... follow up with Jessica and Sarah on this
#    -- intro sign on the S side, introducing what phylogeny is and how it works
#
# steps in the tree (NARROW THIS DOWN TO FIVE KEY MOMENTS):
#  A. All oaks
#  B. OW oaks
#  C. Ilex oaks
#  D. Cerris oaks
#  E1. Ring-cupped oaks
#  E2. Ilex / Cerris oaks
#  F. American oaks
#  G. Red oaks
#  H. Ponticae
#  I. Southern live oaks
#  J. White oaks
#  K. Western white oaks
#  L. White oaks go to eastern North America
#  M. White oaks go back to Eurasia: explain that the position of this node is a bit uncertain, still being studied
#  N. White oaks go down to Mexico [not currently on tree]

sections = list(
  #              "All oaks" = c('gilva', 'alba'),
  #              "B" = c('gilva', 'libani'),
  #              "C" = c('acrodonta', 'coccifera'),
  #              "D" = c('cerris', 'variabilis'),
                "A: Mexican red oaks" = c('delgadoana', 'emoryi'),
                "B: Mexican white oaks" = c('laeta', 'magnoliifolia'),
  #              "Turkey oaks" = c('cerris', 'ilex'),
  #              "Ring-cupped oaks" = c('acuta', 'gilva'),
  #              "American oaks" = c('stellata', 'rubra'),
  #              "G" = c('kelloggii', 'rubra'),
                "C: Ponticae" = c('pontica', 'sadleriana'),
  #              "I" = c('fusiformis', 'minima'),
  #              "White oaks" = c('lobata', 'mongolica'),
  #              "K" = c('lobata', 'garryana'),
  #              "L" = c('alba', 'margaretta'),
                "D: Eurasian white oaks" = c('mongolica', 'petraea')
                )

section.colors = rainbow(length(sections))
background = 'darkgray'
names(section.colors) <- names(sections)

readData <- TRUE
loadTree <- TRUE
relabelTree <- TRUE
relabelNodes <- TRUE
plotTree <- TRUE

labels.geog <- c("Eastern N. America", "Eurasia",
"Mexico, SW U.S., Central America", "Western N. America", "")

if(readData) {
  dat.mor <- read.csv('../DATA/mor.data.geog.garden.v3.csv', as.is = T, row.names = 1)
  codes <- read.csv('../DATA/codes.TSDA.csv', as.is = T, row.names = 1)
  dat.mor$IUCN <- factor(codes[dat.mor$IUCN, 'category'], levels = codes$category)
  dat.mor$Geography <- factor(codes[dat.mor$Geography, 'category'], levels = codes$category)
  dat.mor$Habit <- factor(codes[dat.mor$Habit, 'category'], levels = codes$category)
}

if(loadTree) {
  load('../DATA/tr.uniques.2015-10-15.Rdata')
  if(!exists('tr.full')) {
    tr.full <- root(tr.uniques$sppOnly,
             grep("Castanopsis|Castanea|Lithocarpus", tr.uniques$sppOnly$tip.label))
    #tr.full <- drop.tip(tr.orig, grep("Castanopsis|Castanea|Lithocarpus|Cyclobalanopsis glaucoides", tr$tip.label))
    tr.full <- chronos(tr.full, model = 'discrete')
  }
  tr.full$tip.label <- sapply(strsplit(tr.full$tip.label, " ", fixed = T), function(x) x[2])
#  class(tr.full) <- 'phylo' # added 2018-02-16
  ## added 2018-04-24 below: using New Phyt tree as backbone, welding on eurasian tips
  tr.am <- read.nexus('../DATA/trs.calib.jackknife.4.annotated.tre')
  tr.am <- drop.tip(tr.am, grep('Quercus', tr.am$tip.label, invert=T))
  tr.am <- drop.tip(tr.am, c('Quercus_libani', 'Quercus_trojana'))
  tr.am <- drop.tip(tr.am, c('Quercus_acutissima', 'Quercus_baronii'), trim.internal = FALSE)
  tr.cerris <- extract.clade(tr.full, getMRCA(tr.full, c('acuta', 'coccifera')))
  tcd <- max(node.depth.edgelength(tr.cerris))
  tad <- max(node.depth.edgelength(tr.am)) -
         node.depth.edgelength(tr.am)[which(tr.am$tip.label == 'NA')]
  tr.cerris$edge.length <- tad * (tr.cerris$edge.length / tcd)
  tr <- bind.tree(tr.am, tr.cerris, which(tr.am$tip.label == 'NA'))
  }

if(relabelTree) {
  tr$tip.label <- gsub('Cyclobalanopsis_|Quercus_', '', tr$tip.label)
  tr$tip.label[tr$tip.label == "ithaburense"] <- "ithaburensis"
  tr$tip.label[tr$tip.label == "margarettae"] <- "margaretta"
  tr$tip.label[tr$tip.label == "michauxi"] <- "michauxii"
  tr$tip.label[tr$tip.label == "myrsinaefolia"] <- "myrsinifolia"
  #tr$tip.label[tr$tip.label == "vacciniifolia"] <- "vaccinifolia"
  tr$tip.label[tr$tip.label == "varibilis"] <- "variabilis"
  tr$tip.label[tr$tip.label == 'austrocochinchiensis'] <- 'austrocochinchinensis'
  tr$tip.label[tr$tip.label == 'phillyraeoides'] <- 'phillyreoides'
  tr <- drop.tip(tr, which(!tr$tip.label %in% row.names(dat.mor)[which(dat.mor$use)]))
  tr <- ladderize(tr)
}

if(relabelNodes) {
  tr$node.label <- rep(NA, tr$Nnode + length(tr$tip.label))
  tr.mrca <- mrca(tr)
#  tr.colors <- rep('black', length(tr$edge))
  for(i in names(sections)) {
    tr$node.label[tr.mrca[sections[[i]][1], sections[[i]][2]]] <- i
#    edgesToColor <- which(tr$edge[, 2] %in% getDescendants(tr, tr.mrca[sections[[i]][1], sections[[i]][2]]))
#    tr.colors[edgesToColor] <- section.colors[i]
    }
  tr.colorIndex <- groupClade(tr, node = which(!is.na(tr$node.label))) %>%
    attr('group') %>%
    as.character %>%
    as.numeric
  tr.colors <- section.colors[tr.colorIndex]
  }

## if time -- splice in missing taxa

if(plotTree) {
  tr.plot <- tr
  tr.plot <- groupClade(tr.plot, c(tr.mrca['kelloggii', 'ilicifolia'],
                                   tr.mrca['lobata', 'macranthera'],
                                   tr.mrca['gilva', 'trojana'],
                                   tr.mrca['chrysolepis', 'vacciniifolia'],
                                   tr.mrca['sadleriana', 'pontica'],
                                   tr.mrca['virginiana', 'fusiformis']
                                   )
                                 )
  attr(tr.plot, 'group') <- as.factor(c('All other oaks: white',
                                        'Red oaks: red',
                                        'White oaks: white',
                                        'Eurasian ring-cupped and turkey oaks: black',
                                        'Intermediate oaks: yellow',
                                        'Ponticae: blue',
                                        'Southern live oaks, sect. Virentes: darkgreen')[attr(tr.plot, 'group')])
  attr(tr.plot, 'group')[head(which(attr(tr.plot, 'group') == 'All other oaks: white'),2)] <-
    'Eurasian ring-cupped and turkey oaks: black'
  tr.plot$tip.label <- paste('Q.', tr.plot$tip.label)
  p <- ggtree(tr.plot, layout = 'circular', size = 1.0,
              aes(color = group)
              #color = tr.colors
              #open.angle = 10,
              #color = 'black'
              )
  p <- p + scale_color_manual('Clades',
      values = sapply(strsplit(levels(attr(tr.plot, 'group')), ": ", fixed = T),
                      function(x) x[2])
      )
  p <- p + theme_tree(background)
  p1 <- p + geom_tiplab2(fontface='italic', size = 4,
                         #offset = 2.0)
                         offset = 10,
                         color = 'black')
  dat.geog <- dat.mor[tr$tip.label, c('Geography', 'Habit', 'IUCN')]

  row.names(dat.geog) <- tr.plot$tip.label
  p2 <- gheatmap(p1, dat.geog, width = 0.15,
                font.size = 0,
                #offset = 40.0)
                offset = 0)

  p3 <- p2 + scale_fill_manual('Character states',
                labels = sapply(strsplit(codes$category, "|", fixed = T), function(x) x[2]),
                values = codes$colors
              )
  p3 <- p3 + theme(plot.background = element_rect(fill=background),
                  legend.background = element_rect(fill=background),
                  legend.key = element_rect(fill=background)
                  )

  p4 <- p3 + geom_label(aes(x=branch),
                        label = tr.plot$node.label,
                        size = 3,
                        #fill = section.colors[tr.plot$node.label],
                        fill = 'gray',
                        label.padding = unit(0.18, "lines"),
                        label.r = unit(0.1, "lines"),
                        color = 'black')
  p4b <- p3 + geom_label(aes(x=branch),
                        label = substr(tr.plot$node.label, 1, 1),
                        size = 2,
                        #fill = section.colors[tr.plot$node.label],
                        fill = 'gray',
                        label.padding = unit(0.18, "lines"),
                        label.r = unit(0.1, "lines"),
                        color = 'black')


#   p <- p + theme(legend.position = c(0.15, 0.9),
#                  legend.title = element_text(size = 12),
#                  legend.text = element_text(size = 7),
#                  legend.key.size = unit(0.60, 'cm'),
#                  legend.box.background = element_rect(color = NA)
#                  )
   version <- max(sapply(strsplit(dir(patt = 'tsda'), '[.v]'), function(x) as.numeric(x[3]) + 1))
   pdf(paste("tsdaOaks.v", version, format(Sys.time(), ".%Y-%m-%d.pdf"), sep = ''), 15, 12)
   print(p4)
   print(p4b)
   print(p3)
   dev.off()
   dat.geog.writeVersion <- dat.geog
   for(i in names(dat.geog.writeVersion)) dat.geog.writeVersion[[i]] <- sapply(strsplit(as.character(dat.geog.writeVersion[[i]]), "|", fixed = T), function(x) x[2])
   write.csv(dat.geog.writeVersion, format(Sys.time(), "data.snapshot.for.TSDA.phylogeny.design.%Y-%m-%d.csv"))
}
