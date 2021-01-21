library(PathwayParser)
#library(RUnit)
library(RCyjs)
library(later)
#----------------------------------------------------------------------------------------------------
extractAndWriteNodesAndEdges <- function(id, destinationDirectory="processed")
{
   sbml.filename <- sprintf("%s.sbml", id)
   full.path <- file.path("../incoming", sbml.filename)
   file.size <- file.info(full.path)$size
   #if(file.size > 5000000) {
   #    printf("*** skipping %s, %fMB", id, file.size/1000000)
   #    return()
   #    }
   #printf("--- %s exists? %s", full.path, file.exists(full.path))
   stopifnot(file.exists(full.path))
   pathway <- Pathway$new(full.path)

   reaction.count <- pathway$getReactionCount()
   printf("%s: %d reactions", id, reaction.count)
   tbls <- lapply(seq_len(reaction.count), function(i){
       x <- pathway$processReaction(i, excludeUbiquitousSpecies=TRUE, includeComplexMembers=FALSE)
       })

   tbls.edges <- lapply(tbls, function(el)el$edges)
   tbls.nodes <- lapply(tbls, function(el)el$nodes)
   #checkEquals(length(tbls.edges), reaction.count)
   #checkEquals(length(tbls.nodes), reaction.count)

   tbl.edges <- do.call(rbind, tbls.edges)
   tbl.nodes <- unique(do.call(rbind, tbls.nodes))

   printf("%s: %d nodes, %d edges", id, nrow(tbl.nodes), nrow(tbl.edges))

   #checkEquals(length(unique(c(tbl.edges$source, tbl.edges$target))), nrow(tbl.nodes))


   nodes.file <- file.path(destinationDirectory, sprintf("%s-nodes.tsv", id))
   edges.file <- file.path(destinationDirectory, sprintf("%s-edges.tsv", id))

   write.table(tbl.nodes,
               file=nodes.file,
               sep="\t",
               row.names=FALSE,
               col.names=TRUE,
               quote=FALSE)

   write.table(tbl.edges,
               file=edges.file,
               sep="\t",
               row.names=FALSE,
               col.names=TRUE,
               quote=FALSE)

} # extractAndWriteNodesAndEdges
#----------------------------------------------------------------------------------------------------
run.mtor <- function()
{
   file <- system.file(package="PathwayParser", "extdata", "ReactomePathways-human.tsv")
   tbl.pw <- read.table(file, sep="\t", header=FALSE, as.is=TRUE, nrow=-1, quote="")
   colnames(tbl.pw) <- c("id", "name", "organism")

   ids <- tbl.pw[grep("mtor", tbl.pw$name, ignore.case=TRUE), "id"]
     # 132  R-HSA-9639288                      Amino acids regulate mTORC1 Homo sapiens
     # 770   R-HSA-380972 Energy dependent regulation of mTOR by LKB1-AMPK Homo sapiens
     # 1261  R-HSA-165159                                  MTOR signalling Homo sapiens
     # 2438  R-HSA-166208                       mTORC1-mediated signalling Homo sapiens

   ignore <- lapply(ids, extractAndWriteNodesAndEdges)

} # run.mtor
#----------------------------------------------------------------------------------------------------
run.all <- function()
{
   file <- system.file(package="PathwayParser", "extdata", "ReactomePathways-human.tsv")

   tbl.pw <- read.table(file, sep="\t", header=FALSE, as.is=TRUE, nrow=-1, quote="")
   colnames(tbl.pw) <- c("id", "name", "organism")

   ids <- tbl.pw$id
   ids.big <- rev(c("R-HSA-162582","R-HSA-1430728","R-HSA-168256","R-HSA-1643685","R-HSA-392499",
                    "R-HSA-74160","R-HSA-556833","R-HSA-73857","R-HSA-5663205","R-HSA-168249",
                    "R-HSA-212436","R-HSA-9006934","R-HSA-1280215","R-HSA-597592"))

   print(length(ids.big))
   #ids <- tbl.pw[grep("mtor", tbl.pw$name, ignore.case=TRUE), "id"]

     # 132  R-HSA-9639288                      Amino acids regulate mTORC1 Homo sapiens
     # 770   R-HSA-380972 Energy dependent regulation of mTOR by LKB1-AMPK Homo sapiens
     # 1261  R-HSA-165159                                  MTOR signalling Homo sapiens
     # 2438  R-HSA-166208                       mTORC1-mediated signalling Homo sapiens

   grep("R-HSA-4641262", ids)
   ignore <- lapply(ids.big, function(id) extractAndWriteNodesAndEdges(id, "all"))
   #ignore <- lapply(ids[685:2477], function(id) extractAndWriteNodesAndEdges(id, "all"))

} # run.all
#----------------------------------------------------------------------------------------------------
view <- function()
{
   if(!exists("rcy")){
     rcy <<- RCyjs()
     setBrowserWindowTitle(rcy, "mTOR-related")
     }


   reactions <- c("R-HSA-392499", "R-HSA-5368286", "R-HSA-5368287", "R-HSA-5389840",
                  "R-HSA-5419276", "R-HSA-72766")

   subset(tbl.pw, id %in% reactions)$name

   reaction <- reactions[1]  # edges: 3854   nodes: 2931        "Metabolism of proteins"
   reaction <- reactions[2]  # edges: 21   nodes: 20            "Mitochondrial translation"
   reaction <- reactions[3]  # edges: 76   nodes: 58            "Mitochondrial translation elongation"
   reaction <- reactions[4]  # edges: 30   nodes: 26            "Mitochondrial translation initiation"
   reaction <- reactions[5]  # edges: 25   nodes: 22            "Mitochondrial translation termination"
   reaction <- reactions[6]  # edges: 525   nodes: 387          "Translation"

   deleteGraph(rcy)

   for(reaction in reactions[2:5]){
      edges.filename <- sprintf("all-19jan2021/%s-edges.tsv", reaction)
      nodes.filename <- sprintf("all-19jan2021/%s-nodes.tsv", reaction)
      file.exists(edges.filename)
      file.exists(nodes.filename)
      tbl.edges <- read.table(edges.filename, sep="\t", header=TRUE, as.is=TRUE, quote="")
      tbl.nodes <- read.table(nodes.filename, sep="\t", header=TRUE, as.is=TRUE, quote="")
      printf("edges: %d   nodes: %d", nrow(tbl.edges), nrow(tbl.nodes))
      g.json <- toJSON(dataFramesToJSON(tbl.edges, tbl.nodes))
      addGraph(rcy, g.json)
      } # for reaction

   later(function(){
      setBrowserWindowTitle(rcy, sprintf("%s", reaction))
      loadStyleFile(rcy, "style.js")
      layout(rcy, "cola") # "cose-bilkent")
      fit(rcy)
      }, 2)


} # view
#----------------------------------------------------------------------------------------------------
