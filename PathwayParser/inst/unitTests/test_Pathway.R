library(PathwayParser)
library(RUnit)
library(RCyjs)
library(EnsDb.Hsapiens.v79)
library(later)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_getReactionNames()
   test_getPathways()
   test_getMolecularSpeciesMap()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   sbml.filename <- system.file(package="PathwayParser", "extdata", "R-HSA-165159.sbml")
   checkTrue(file.exists(sbml.filename))
   pathway <- Pathway$new(sbml.filename)
   checkTrue("Pathway" %in% is(pathway))
   checkEquals(pathway$getName(), "MTOR signalling")
   checkEquals(pathway$getReactionCount(), 29)

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_getReactionNames <- function()
{
   message(sprintf("--- test_getReactionNames"))
   sbml.filename <- system.file(package="PathwayParser", "extdata", "R-HSA-165159.sbml")
   pathway <- Pathway$new(sbml.filename)
   names <- pathway$getReactionNames()
   checkEquals(length(names), pathway$getReactionCount())
       # spot check:
   checkTrue(all(c("Ragulator binds Rag dimers",
                   "Formation of active mTORC1 complex",
                   "Rag dimer formation",
                   "AKT1S1 (PRAS40) binds mTORC1") %in% names))


} # test_getReactionNames
#------------------------------------------------------------------------------------------------------------------------
test_getPathways <- function()
{
   message(sprintf("--- test_getPathways"))
   sbml.filename <- system.file(package="PathwayParser", "extdata", "R-HSA-165159.sbml")
   pathway <- Pathway$new(sbml.filename)
   tbl.pathways <- pathway$getPathwayNameMap()
   checkEquals(ncol(tbl.pathways), 3)
   checkEquals(colnames(tbl.pathways), c("id", "name", "organism"))
   checkTrue(nrow(tbl.pathways) > 2000)
   checkTrue(nrow(tbl.pathways) < 3000)
   checkTrue(all(tbl.pathways$organism == "Homo sapiens"))

} # test_getPathways
#------------------------------------------------------------------------------------------------------------------------
test_getMolecularSpeciesMap <- function()
{
   message(sprintf("--- test_getMolecularSpeciesMap"))
   sbml.filename <- system.file(package="PathwayParser", "extdata", "R-HSA-165159.sbml")
   pathway <- Pathway$new(sbml.filename)
   tbl.species <- pathway$getMolecularSpeciesMap()
   checkEquals(dim(tbl.species), c(66, 5))
   checkEquals(colnames(tbl.species), c("id", "name", "type", "compartment", "complex.members"))
   checkEquals(unlist(lapply(tbl.species, class), use.names=FALSE),
               c("character", "character", "character", "character", "list"))
   checkEquals(tbl.species$id[1], "species_5653921")
   checkEquals(tbl.species$name[1], "Ragulator")
   checkEquals(tbl.species$type[1], "complex")
   checkEquals(tbl.species$compartment[1], "lysosomal membrane")
   checkEquals(tbl.species$complex.members[1][[1]], c("uniprotkb:O43504", "uniprotkb:Q6IAA8",
                                                      "uniprotkb:Q9Y2Q5", "uniprotkb:Q9UHA4",
                                                      "uniprotkb:Q0VGL1"))

} # test_getMolecularSpeciesMap
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()