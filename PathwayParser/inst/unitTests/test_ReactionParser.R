library(PathwayParser)
library(RUnit)
library(RCyjs)
library(EnsDb.Hsapiens.v79)
library(later)
#------------------------------------------------------------------------------------------------------------------------
# 2nd reaction to test, #3
#   <reaction compartment="compartment_17938"  id="reaction_165718"
#      name="mTORC1 phosphorylation of RPS6KB1 (S6K)" reversible="false">
#------------------------------------------------------------------------------------------------------------------------
sbml.filename <- system.file(package="PathwayParser", "extdata", "R-HSA-165159.sbml")
if(!exists("doc")){
   stopifnot(file.exists(sbml.filename))
   text <- paste(readLines(sbml.filename), collapse="\n")
   checkTrue(nchar(text) > 300000)   # 315493
   doc <- read_xml(text)
   xml_ns_strip(doc)
   }
pathway <- Pathway$new(sbml.filename)
#------------------------------------------------------------------------------------------------------------------------
# note that the whole docuent is returned, but some mysterious internal state records
# the selection made here.  good to check this by calling xml_path, which should
# return /sbml/model/listOfReactions/reaction[i]
getReactionForTesting <- function(doc, i)
{
    xml_find_all(doc, sprintf("//reaction[%d]", i))

} # getReactionForTesting
#------------------------------------------------------------------------------------------------------------------------
if(!exists("reaction"))
   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
if(!exists("parser"))
   parser <- ReactionParser$new(doc, reaction, pathway)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
  # test_pathwayNameMap()
   test_getReactants()
   test_getProducts()
   test_getModifiers()
   #test_molecularSpeciesMap()
   test_eliminateUbiquitiousSpecies()
   test_toEdgeAndNodeTables()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))
   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)
   checkEquals(xml_path(reaction), "/sbml/model/listOfReactions/reaction[5]")

   checkEquals(parser$getXPath(), "/sbml/model/listOfReactions/reaction[5]")
   checkEquals(parser$getID(), "reaction_9679044")
   checkEquals(parser$getName(), "FKBP1A binds sirolimus")
   checkEquals(parser$getCompartment(), "compartment_70101")
   notes <- parser$getNotes()
   checkTrue(nchar(notes) > 1100)
   checkEquals(substring(notes, 1, 33), "Sirolimus is a macrolide compound")

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
# test_pathwayNameMap <- function()
# {
#    parser <- ReactionParser$new(doc, reaction, pathway)
#    tbl.pathwayNames <- parser$getPathwayNameMap()
#    checkEquals(colnames(tbl.pathwayNames), c("id", "name", "organism"))
#    checkTrue(nrow(tbl.pathwayNames) > 2400)
#    checkEquals(subset(tbl.pathwayNames, name=="Glycolysis")$id, "R-HSA-70171")
#
#    f <- system.file(package="PathwayParser", "extdata", "R-HSA-70171.sbml")
#    checkTrue(file.exists(f))
#    text <- paste(readLines(f), collapse="\n")
#    doc2 <- read_xml(text)
#    xml_ns_strip(doc2)
#
#    reaction <- getReactionForTesting(doc2, 1)
#    parser2 <- ReactionParser$new(doc2, reaction, pathway)
#    parser2$toEdgeAndNodeTables(T,F)
#    parser2$getName()
#    length(xml_find_all(doc2, "//reaction")) # 24 reactions in glycolysis
#    xml_text(xml_find_all(doc2, "//reaction/@name")) # 24 reactions in glycolysis
#
#      # [1] "D-fructose 6-phosphate + ATP => D-fructose 1,6-bisphosphate + ADP"
#      # [2] "ADPGK:Mg2+ phosphorylates Glc to G6P"
#      # [3] "nucleoplasmic GCK1:GKRP complex => glucokinase (GCK1) + glucokinase regulatory protein (GKRP)"
#      # [4] "glucokinase (GCK1) + glucokinase regulatory protein (GKRP) <=> GCK1:GKRP complex"
#      # [5] "glucokinase [nucleoplasm] => glucokinase [cytosol]"
#      # [6] "cytosolic GCK1:GKRP complex <=> glucokinase (GCK1) + glucokinase regulatory protein (GKRP)"
#      # [7] "NPC transports GCK1:GKRP from cytosol to nucleoplasm"
#      # [8] "alpha-D-glucose 6-phosphate <=> D-fructose 6-phosphate"
#      # [9] "3-Phospho-D-glycerate <=> 2-Phospho-D-glycerate"
#      #[10] "D-fructose 1,6-bisphosphate <=> dihydroxyacetone phosphate + D-glyceraldehyde 3-phosphate"
#      #[11] "1,3-bisphospho-D-glycerate + ADP <=> 3-phospho-D-glycerate + ATP"
#      #[12] "2-Phospho-D-glycerate <=> Phosphoenolpyruvate + H2O"
#      #[13] "D-glyceraldehyde 3-phosphate + orthophosphate + NAD+ <=> 1,3-bisphospho-D-glycerate + NADH + H+"
#      #[14] "GNPDA1,2 hexamers deaminate GlcN6P to Fru(6)P"
#      #[15] "HK1,2,3,GCK phosphorylate Glc to form G6P"
#      #[16] "dihydroxyacetone phosphate <=> D-glyceraldehyde 3-phosphate"
#      #[17] "phosphoenolpyruvate + ADP => pyruvate + ATP"
#      #[18] "BPGM dimer isomerises 1,3BPG to 2,3BPG"
#      #[19] "PGM2L1:Mg2+ phosphorylates G6P to G1,6BP"
#      #[20] "PGP:Mg2+ dimer hydrolyses 3PG to glycerol"
#      #[21] "Dephosphorylation of phosphoPFKFB1 by PP2A complex"
#      #[22] "Fructose 2,6-bisphosphate is hydrolyzed to form fructose-6-phosphate and orthophosphate"
#      #[23] "D-fructose 6-phosphate + ATP => D-fructose 2,6-bisphosphate + ADP"
#      #[24] "Phosphorylation of PF2K-Pase by PKA catalytic subunit"
#
#
# } # test_pathwayNameMap
#------------------------------------------------------------------------------------------------------------------------
test_getReactants <- function()
{
   message(sprintf("--- test_getReactants"))

   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)

   checkEquals(parser$getReactantCount(), 2)
   checkEquals(sort(parser$getReactants()), sort(c("species_2026007", "species_9678687")))

} # test_getReactants
#------------------------------------------------------------------------------------------------------------------------
test_getProducts <- function()
{
   message(sprintf("--- test_getProducts"))

   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)

   checkEquals(parser$getProductCount(), 1)
   checkEquals(sort(parser$getProducts()), "species_9679098")

} # test_getProducts
#------------------------------------------------------------------------------------------------------------------------
# reaction 2 in mTOR signaling hsa a modifier
#   "Phosphorylation and activation of eIF4B by activated S6K1"
# make sure we can find it here
test_getModifiers <- function()
{
   message(sprintf("--- test_getModifiers"))

   reaction.2 <- getReactionForTesting(doc, 2)
   parser.2 <- ReactionParser$new(doc, reaction.2, pathway)
   checkEquals(parser.2$getModifierCount(), 1)
   checkEquals(sort(parser.2$getModifiers()), "species_165714")

     # get the name - though it is an awkward one, not the simple "S6K1" we'd like
     # wikipedia:
     #   Ribosomal protein S6 kinase beta-1 (S6K1),
     #   also known as p70S6 kinase (p70S6K, p70-S6K)
     #   The phosphorylation of p70S6K at threonine 389 has been used as a hallmark of activation
     #   by mTOR and correlated with autophagy inhibition in various situations. However, several
     # recent studies suggest that the activity of p70S6K plays a more positive role in the increase
     # of autophagy.

   tbl.map <- pathway$getMolecularSpeciesMap()
   checkEquals(subset(tbl.map, id == "species_165714")$name, "p-S371,T389-RPS6KB1")
   checkEquals(subset(tbl.map, id == "species_165714")$type, "molecule")
   checkEquals(subset(tbl.map, id == "species_165714")$compartment, "cytosol")
   checkTrue(is.null(subset(tbl.map, id == "species_165714")$members))

} # test_getProducts
#------------------------------------------------------------------------------------------------------------------------
# test_molecularSpeciesMap <- function()
# {
#    message(sprintf("--- test_molecularSpeciesMap"))
#
#   x <- parser$getMolecularSpeciesMap()
#    checkTrue(is.list(x))
#    checkEquals(length(x), 66)
#    moleculeTypes <- unlist(lapply(x, "[", "moleculeType"))
#    counts <- as.list(table(moleculeTypes))
#
#    checkEquals(counts$complex, 30)
#    checkEquals(counts$molecule, 36)
#
# } # test_getProducts
#------------------------------------------------------------------------------------------------------------------------
test_toEdgeAndNodeTables <- function()
{
   message(sprintf("--- test_toEdgeAndNodeTables"))
      # user reaction 2 in R-HSA-165159.sbml: 1 each of reactant, product, modifier
   reaction.2 <- getReactionForTesting(doc, 2)
   parser.tmp <- ReactionParser$new(doc, reaction.2, pathway)
   checkEquals(parser.tmp$getReactantCount(), 2)
   checkEquals(parser.tmp$getProductCount(), 2)
   checkEquals(parser.tmp$getModifierCount(), 1)
   checkEquals(sort(parser.tmp$getModifiers()), "species_165714")

   x <- parser.tmp$toEdgeAndNodeTables(includeComplexMembers=FALSE)
   checkEquals(sort(names(x)), c("edges", "nodes"))

   checkEquals(dim(x$nodes), c(4, 4))
   checkTrue("species_165714" %in% x$nodes$id)
   checkTrue("p-S371,T389-RPS6KB1" %in% x$nodes$label)

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x$edges), c(3, 3))
   checkEquals(x$edges[3, "source"], "species_165714")
   checkEquals(x$edges[3, "target"], "reaction_165777")
   checkEquals(x$edges[3, "interaction"], "modifies")

      # now get all species, including water, atp, adp if present
   x <- parser.tmp$toEdgeAndNodeTables(includeComplexMembers=FALSE, excludeUbiquitousSpecies=FALSE)
   checkEquals(sort(names(x)), c("edges", "nodes"))

   checkEquals(dim(x$nodes), c(6, 4))
   checkTrue(all(c("ATP", "ADP") %in% x$nodes$label))

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x$edges), c(5, 3))

} # test_toEdgeAndNodeTables
#------------------------------------------------------------------------------------------------------------------------
renderReaction <- function()
{
   x <- parser$toEdgeAndNodeTables()
   g.json <- toJSON(dataFramesToJSON(x$edges, x$nodes))
   deleteGraph(rcy)
   addGraph(rcy, g.json)
   loadStyleFile(rcy, "style.js")
   layout(rcy, "cose")
   fit(rcy)

} # renderReaction
#------------------------------------------------------------------------------------------------------------------------
test_reaction_1 <- function()
{
   message(sprintf("--- test_reaction_1"))
   reaction <- getReactionForTesting(doc, 1) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$toEdgeAndNodeTables()

   g.json <- toJSON(dataFramesToJSON(x$edges, x$nodes))
   deleteGraph(rcy)
   addGraph(rcy, g.json)
   loadStyleFile(rcy, "style.js")
   layout(rcy, "cose")
   fit(rcy)

} # test_reaction_1
#------------------------------------------------------------------------------------------------------------------------
test_eliminateUbiquitiousSpecies <- function()
{
   message(sprintf("--- test_eliminateUbiquitousSpecies"))

   reaction <- getReactionForTesting(doc, 3) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$toEdgeAndNodeTables(excludeUbiquitousSpecies=FALSE, includeComplexMembers=FALSE)
   checkEquals(nrow(x$edges), 8)
   checkEquals(nrow(x$nodes), 9)

   x <- parser$toEdgeAndNodeTables(excludeUbiquitousSpecies=TRUE, includeComplexMembers=FALSE)
   checkEquals(nrow(x$edges), 6)
   checkEquals(nrow(x$nodes), 7)

   x <- parser$toEdgeAndNodeTables(excludeUbiquitousSpecies=FALSE, includeComplexMembers=FALSE)
   checkEquals(nrow(x$edges), 8)
   checkEquals(nrow(x$nodes), 9)

} # test_eliminateUbiquitousSpecies
#------------------------------------------------------------------------------------------------------------------------
displayReaction <- function(i, exclude=TRUE, deleteExistingGraph=TRUE, includeComplexMembers=FALSE)
{
   if(!exists("rcy")){
      rcy <<- RCyjs()
      setBrowserWindowTitle(rcy, "ReactionParser")
      }

   reaction <- getReactionForTesting(doc, i)
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$toEdgeAndNodeTables(excludeUbiquitousSpecies=TRUE, includeComplexMembers)

   g.json <- toJSON(dataFramesToJSON(x$edges, x$nodes))

   if(deleteExistingGraph)
      deleteGraph(rcy)

   addGraph(rcy, g.json)
   later(function(){
      setBrowserWindowTitle(rcy, sprintf("%d: %s", i, parser$getName()))
      loadStyleFile(rcy, "style.js")
      layout(rcy, "cose-bilkent")
      fit(rcy)
      }, 2)

} # displayReaction
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
