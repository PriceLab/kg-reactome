library(PathwayParser)
library(RUnit)
library(RCyjs)
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
   test_getCounts()
   test_getReactants()
   test_getProducts()
   test_getModifiers()
   #test_getComplexes()
   test_assignNodeType()
   test_assignNodeName()
   # test_eliminateUbiquitiousSpecies()
   test_toEdgeAndNodeTables()
   #test_toEdgeAndNodeTables_withComplexes()

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
test_assignNodeType <- function()
{
   message(sprintf("--- test_assignNodeType"))
   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)
   species <- unique(c(parser$getReactants(), parser$getProducts(), parser$getModifiers(),
                unlist(parser$getComplexes(), use.names=FALSE),
                "species_9678687"))
   x <- lapply(species, parser$assignNodeType)
   names(x) <- species

     # TODO - more explict test needed, should be restored.
   checkEquals(sort(as.character(x)), c("complex","drug", "ligand", "protein", "protein"))
   #checkEquals(x, list(species_9678687="drug",        # sirolimus [cytosol]
   #                    species_2026007="protein",     # FKB1A [cytosol]
   #                    species_9679098="complex",     # FKBP1A:sirolimus [cytosol]
   #                    `uniprotkb:P62942`="protein",  # FKBP1A
   #                    `ligandId:6031`="ligand",      # another identifier for sirolimus (rapamycin)
   #                    species_9678687="drug"))       # added explicitly, an earlier test, redundant?

} # test_assignNodeType
#------------------------------------------------------------------------------------------------------------------------
test_assignNodeName <- function()
{
   message(sprintf("--- test_assignNodeName"))
   reaction <- getReactionForTesting(doc, 5) #    # "FKBP1A binds sirolimus"
   parser <- ReactionParser$new(doc, reaction, pathway)
   species <- unique(c(parser$getReactants(), parser$getProducts(), parser$getModifiers(),
                unlist(parser$getComplexes(), use.names=FALSE)))
#                "species_9678687")
   x <- lapply(species, parser$assignNodeName)
   names(x) <- species

   checkEquals(x$species_9678687, "rapamycin")
   checkEquals(x$species_2026007, "FKBP1A")
   checkEquals(x$species_9679098, "FKBP1A:sirolimus")
   checkEquals(x$`uniprotkb:P62942`, "FKBP1A")
   checkEquals(x$`ligandId:6031`, "rapamycin")

   #checkEquals(x, list(species_9678687="rapamycin",
   #                    species_2026007="FKBP1A",
   #                    species_9679098="FKBP1A:sirolimus",
   #                    `uniprotkb:P62942`="FKBP1A",
   #                    `ligandId:6031`="rapamycin",
   #                    species_9678687="rapamycin"))

} # test_assignNodeName
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

      # a member of the ribosomal S6 kinase family of serine/threonine
      # kinases. The encoded protein responds to mTOR (mammalian
      # target of rapamycin) signaling to promote protein synthesis,
      # cell growth, and cell proliferation.
      # here phosphorylated on serine 371, threonine 389
      # phosphorylating RPS6 in reaction R-HSA-165726

   x <- as.list(subset(tbl.map, id == "species_165714"))
   checkEquals(x$id, "species_165714")
   checkEquals(x$name, "p-S371,T389-RPS6KB1")
   checkEquals(x$type, "protein")
   checkEquals(x$compartment, "cytosol")
   checkEquals(x$uniprot.id, "uniprotkb:P23443")
   checkTrue(is.na(x$chebi.id))
   checkEquals(x$complex.members[[1]], "uniprotkb:P23443")   # does this make sense?  TODO

      # ATP
   x <- as.list(subset(tbl.map, id == "species_113592"))
   with(x,
        checkEquals(id, "species_113592"),
        checkEquals(name, "ATP"),
        checkEquals(type,"smallMolecule"),
        checkEquals(compartment, "cytosol"),
        checkTrue(is.na(uniprot.id)),
        checkEquals(chebi.id, "CHEBI:30616"),
        checkEquals(complex.members, list(NULL))
        )

} # test_getModifiers
#------------------------------------------------------------------------------------------------------------------------
test_getComplexes <- function()
{
   message(sprintf("--- test_getComplexes"))

      #-----------------------------------------------------------------
      # use reaction 5 in R-HSA-165159.sbml:  1 complex of two elements
      #-----------------------------------------------------------------

   reaction.5 <- getReactionForTesting(doc, 5)
   parser <- ReactionParser$new(doc, reaction.5, pathway)
   complex.list <- parser$getComplexes()
   checkEquals(complex.list, list(species_9679098=c("uniprotkb:P62942", "ligandId:6031")))

      #-----------------------------------------------------------------
      # use reaction 2 in R-HSA-165159.sbml:  no complexes
      #-----------------------------------------------------------------

   reaction.2 <- getReactionForTesting(doc, 2)
   parser <- ReactionParser$new(doc, reaction.2, pathway)
   complex.list <- parser$getComplexes()
   checkEquals(length(complex.list), 0)

      #-----------------------------------------------------------------
      # use reaction 1 in R-HSA-165159.sbml: 3 complexes
      #-----------------------------------------------------------------

   reaction.1 <- getReactionForTesting(doc, 1)
   parser <- ReactionParser$new(doc, reaction.1, pathway)
   complex.list <- parser$getComplexes()
   checkEquals(complex.list,
               list(species_5653945=c("uniprotkb:Q9HB90", "ChEBI:17552", "uniprotkb:Q9NQL2",
                                       "ChEBI:15996", "uniprotkb:Q7L523", "uniprotkb:Q5VZM2"),
                    species_5653921=c("uniprotkb:O43504", "uniprotkb:Q6IAA8", "uniprotkb:Q9Y2Q5",
                                      "uniprotkb:Q9UHA4", "uniprotkb:Q0VGL1"),
                    species_5653979=c("uniprotkb:O43504", "uniprotkb:Q6IAA8", "uniprotkb:Q9Y2Q5",
                                      "uniprotkb:Q9UHA4", "uniprotkb:Q0VGL1", "uniprotkb:Q9HB90",
                                      "ChEBI:17552", "uniprotkb:Q9NQL2", "ChEBI:15996",
                                      "uniprotkb:Q7L523", "uniprotkb:Q5VZM2")))


} # test_getComplexes
#------------------------------------------------------------------------------------------------------------------------
test_getCounts <- function()
{
   message(sprintf("--- test_getCounts"))

   reaction <- getReactionForTesting(doc, 1)
   parser <- ReactionParser$new(doc, reaction, pathway)
   checkEquals(parser$getReactantCount(), 2)
   checkEquals(parser$getProductCount(), 1)
   checkEquals(parser$getModifierCount(), 0)
   checkEquals(parser$getComplexCount(), 3)

   reaction <- getReactionForTesting(doc, 11)
   parser <- ReactionParser$new(doc, reaction, pathway)
   checkEquals(parser$getReactantCount(), 5)
   checkEquals(parser$getProductCount(), 1)
   checkEquals(parser$getModifierCount(), 0)
   checkEquals(parser$getComplexCount(), 6)

   reaction <- getReactionForTesting(doc, 29)
   parser <- ReactionParser$new(doc, reaction, pathway)
   checkEquals(parser$getReactantCount(), 2)
   checkEquals(parser$getProductCount(), 2)
   checkEquals(parser$getModifierCount(), 1)
   checkEquals(parser$getComplexCount(), 4)

   reaction.count <- length(xml_find_all(doc, "//reaction"))

   for(i in seq_len(reaction.count)){
     reaction <- getReactionForTesting(doc, i)
     parser <- ReactionParser$new(doc, reaction, pathway)
     checkTrue(parser$getReactantCount() > 0)
     checkTrue(parser$getProductCount() > 0)
     checkTrue(parser$getModifierCount() >= 0)
     checkTrue(parser$getComplexCount() >= 0)
     #printf("%d: %d %d %d %d", i, parser$getReactantCount(),
     #       parser$getProductCount(),
     #       parser$getModifierCount(),
     #       parser$getComplexCount())
     } # for i

      #-------------------------------------------------------
      # counts for all categories, reaction 1
      # "Ragulator binds Rag dimers"
      #-------------------------------------------------------
   reaction <- getReactionForTesting(doc, 1)
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$getCounts()
   checkEquals(x, list(reactants=2, products=1, modifiers=0, complexes=3))

      #-------------------------------------------------------
      # counts for all categories, reaction 8
      # "Phosphorylation of 4E-BP1 by activated mTORC1"
      #-------------------------------------------------------
   reaction <- getReactionForTesting(doc, 8)
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$getCounts()
   checkEquals(x, list(reactants=2, products=2, modifiers=4, complexes=4))

} # test_getCounts
#------------------------------------------------------------------------------------------------------------------------
# reaction 2 in R-HSA-165159.sbml: no complexes
# reaction 5 in R-HSA-165159.sbml: 1 complex of two elements
# reaction 1 in R-HSA-165159.sbml: 3 complexes
test_toEdgeAndNodeTables <- function()
{
   message(sprintf("--- test_toEdgeAndNodeTables"))

      #-----------------------------------------------------------------------------
      # user reaction 2 in R-HSA-165159.sbml: 1 each of reactant, product, modifier
      # but note: there are no complexes in this reaction.
      #-----------------------------------------------------------------------------

   reaction.2 <- getReactionForTesting(doc, 2)
   parser.tmp <- ReactionParser$new(doc, reaction.2, pathway)
   checkEquals(parser.tmp$getReactantCount(), 2)
   checkEquals(parser.tmp$getProductCount(), 2)
   checkEquals(parser.tmp$getModifierCount(), 1)
   checkEquals(sort(parser.tmp$getModifiers()), "species_165714")
   checkEquals(length(parser.tmp$getComplexes()), 4)

     #--------------------------------------------------------------------------
     # do not request inclusion of members of any complex
     #--------------------------------------------------------------------------

   x <- parser.tmp$toEdgeAndNodeTables(includeComplexMembers=FALSE)
   checkEquals(sort(names(x)), c("edges", "nodes"))

     # do some node checks
   checkEquals(dim(x$nodes), c(4, 4))
   checkTrue("species_165714" %in% x$nodes$id)
   checkTrue("p-S371,T389-RPS6KB1" %in% x$nodes$label)
     # check the node types
   checkEquals(x$nodes$id, c("species_72589", "reaction_165777", "species_165714", "species_165773"))
   checkEquals(x$nodes$type, c("protein", "reaction", "protein", "protein"))
   checkEquals(x$nodes$parent, c("", "", "", ""))  # includeComplexMebers=FALSE

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x$edges), c(3, 4))
   checkEquals(x$edges[3, "source"], "species_165714")
   checkEquals(x$edges[3, "target"], "reaction_165777")
   checkEquals(x$edges[3, "reaction"], "Phosphorylation and activation of eIF4B by activated S6K1")

     #--------------------------------------------------------------------------
     # now request inclusion of members of any complex - of which there are none
     #--------------------------------------------------------------------------

   x2 <- parser.tmp$toEdgeAndNodeTables(includeComplexMembers=TRUE)
   checkEquals(x, x2)

      # now get all species, including water, atp, adp if present
   x3 <- parser.tmp$toEdgeAndNodeTables(includeComplexMembers=FALSE, excludeUbiquitousSpecies=FALSE)
   checkEquals(sort(names(x3)), c("edges", "nodes"))

   checkEquals(dim(x3$nodes), c(6, 4))
   checkTrue(all(c("ATP", "ADP") %in% x3$nodes$label))

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x3$edges), c(5, 4))

} # test_toEdgeAndNodeTables
#----------------------------------------------------------------------------------------------------
# reaction 21 has 3 reactant proteins, but two have been interpreted as complexes, since
# they have two alternate members.
# this seems like a funky use of sbml, but let's work around it, construing the multiple
# alternate proteins as reactants to the reaction into which STRAD feeds
# this adjustment will happen in the toEdgeAndNodeTables method
test_toEdgeAndNodeTables_STRAD <- function()
{
   message(sprintf("--- test_toEdgeAndNodeTables"))

   reaction <- getReactionForTesting(doc, 21)
   parser <- ReactionParser$new(doc, reaction, pathway)
   parser$getReactantCount()
   parser$getProductCount()
   parser$getModifierCount()
   parser$getModifiers()
   parser$getComplexes()

     #--------------------------------------------------------------------------
     # do not request inclusion of members of any complex
     #--------------------------------------------------------------------------

   x <- parser$toEdgeAndNodeTables(includeComplexMembers=FALSE)
   checkEquals(sort(names(x)), c("edges", "nodes"))

     # do some node checks
   checkEquals(dim(x$nodes), c(4, 4))
   checkTrue("species_165714" %in% x$nodes$id)
   checkTrue("p-S371,T389-RPS6KB1" %in% x$nodes$label)
     # check the node types
   checkEquals(x$nodes$id, c("species_72589", "reaction_165777", "species_165714", "species_165773"))
   checkEquals(x$nodes$type, c("protein", "reaction", "protein", "protein"))
   checkEquals(x$nodes$parent, c("", "", "", ""))  # includeComplexMebers=FALSE

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x$edges), c(3, 4))
   checkEquals(x$edges[3, "source"], "species_165714")
   checkEquals(x$edges[3, "target"], "reaction_165777")
   checkEquals(x$edges[3, "reaction"], "Phosphorylation and activation of eIF4B by activated S6K1")

     #--------------------------------------------------------------------------
     # now request inclusion of members of any complex - of which there are none
     #--------------------------------------------------------------------------

   x2 <- parser$toEdgeAndNodeTables(includeComplexMembers=TRUE)
   checkEquals(x, x2)

      # now get all species, including water, atp, adp if present
   x3 <- parser$toEdgeAndNodeTables(includeComplexMembers=FALSE, excludeUbiquitousSpecies=FALSE)
   checkEquals(sort(names(x3)), c("edges", "nodes"))

   checkEquals(dim(x3$nodes), c(6, 4))
   checkTrue(all(c("ATP", "ADP") %in% x3$nodes$label))

     # keep in mind that excludeUbiquitousSpecies is default TRUE
   checkEquals(dim(x3$edges), c(5, 4))

} # test_toEdgeAndNodeTables_STRAD
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
# a reaction with all four kinds of nodes
# "Phosphorylation of 4E-BP1 by activated mTORC1"
#   reactants:  2
#   products:   2
#   modifiers:  4
#   complexes:  2
test_reaction_8 <- function(display=FALSE)
{
   message(sprintf("--- test_reaction_8"))

   reaction <- getReactionForTesting(doc, 8) #
   parser <- ReactionParser$new(doc, reaction, pathway)
   x <- parser$toEdgeAndNodeTables(includeComplexMembers=FALSE)

   checkEquals(lapply(x, dim), list(edges=c(6,4), nodes=c(7,4)))

   x <- parser$toEdgeAndNodeTables(includeComplexMembers=TRUE)
   checkEquals(lapply(x, dim), list(edges=c(10,4), nodes=c(9,4)))


} # test_reaction_8
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
# intersecting genehancer and human accelerated regions turns up this mitochondrial protein
# make sure we extract it properly
test_oax1l.Q15070.inR_HSA_72766.sbml <- function()
{
   message(sprintf("--- test_oax1l.Q15070.inR_HSA_72766.sbml"))
   sbml.filename <- system.file("PathwayParser", "extdata", "R-HSA-72766.sbml")
   text <- paste(readLines(sbml.filename), collapse="\n")
   checkTrue(nchar(text) > 300000)   # 315493
   doc <- read_xml(text)
   xml_ns_strip(doc)


} # test_oax1l.Q15070.inR_HSA_72766.sbml
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

   invisible(x)

} # displayReaction
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
