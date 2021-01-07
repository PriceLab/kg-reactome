# Pathway.R
#
#' import R6
#' import xml2
#' import data.table
#'
#' @title Pathway
#' @description an R6 class which parses and represents a Reactome Pathway
#' @name Pathway
#'
#' @export

Pathway = R6Class("Pathway",

    private = list(
       sbml.filename = NULL,
       doc = NULL,
       reactions = NULL,
       tbl.molecularSpecies = NULL,
       tbl.pathwayNameMap = NULL,

       readPathwayIdMap = function(){
          file <- system.file(package="PathwayParser", "extdata", "ReactomePathways-human.tsv")
          tbl.pathwayNameMap <- read.table(file, sep="\t", header=FALSE, as.is=TRUE, nrow=-1, quote="")
          colnames(tbl.pathwayNameMap) <- c("id", "name", "organism")
          private$tbl.pathwayNameMap <- tbl.pathwayNameMap
          },

        createMolecularSpeciesTable = function(){
           tmp <- list()
           all.species <- xml_find_all(private$doc, "..//listOfSpecies/species")
           for(species in all.species){
              id <- xml_text(xml_find_all(species, ".//@id"))
              name.raw <- xml_text(xml_find_all(species, ".//@name"))
              tokens <- strsplit(name.raw, " [", fixed=TRUE)[[1]]
              name <- tokens[1]
              moleculeType <- "molecule"
              compartment <- sub("]", "", tokens[2], fixed=TRUE)
              members <- xml_find_all(species, ".//bqbiol:hasPart//rdf:li/@rdf:resource")
              if(length(members) == 0)
                members <- c()
              if(length(members) > 0){
                moleculeType <- "complex"
                members <- xml_text(members)
                members <- sub("http://purl.uniprot.org/uniprot/", "uniprotkb:", members, fixed=TRUE)
                members <- sub("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:", "ChEBI:", members, fixed=TRUE)
                members <- sub("http://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=",
                               "ligandId:", members, fixed=TRUE)
                }
              new.entry <- list(name=name, moleculeType=moleculeType, compartment=compartment, members=members)
              tmp[[id]] <- new.entry
              } # for species
           extract <- function(i){
               x <- tmp[[i]]
               complex.members <- x$members
               data.table(id=names(tmp)[i], name=x$name, type=x$moleculeType,
                          compartment=x$compartment, complex.members=list(complex.members))
               }
           tbls <- lapply(seq_len(length(tmp)), extract)
           private$tbl.molecularSpecies <- do.call(rbind, tbls)
           } # createMolecularSpeciesTable

       ),

    public = list(

        #' @description
        #' Create a new Pathway object
        #' @param sbml.filename character string, path to a valid sbml 3.1 document
        #' @param xml.node an XMLInternalNode
        #' @return A new `ReactionParser` object.

       initialize = function(sbml.filename){
          stopifnot(file.exists(sbml.filename))
          text <- paste(readLines(sbml.filename), collapse="\n")
          private$doc <- read_xml(text)
          xml_ns_strip(private$doc)
          private$createMolecularSpeciesTable()
          private$readPathwayIdMap()
          },

        #' @description
        #' every Reactome pathway has a name, e.g., R-HSA-165159 is "MTOR signalling"
        #' @return a character string, the pathway name

       getName = function(){
          xml_text(xml_find_all(private$doc, "//model/@name"))
          },

        #' @description
        #' every Reactome pathway has multiple reactions
        #' @return numeric count
       getReactionCount = function(){
          length(xml_find_all(private$doc, "//reaction"))

          },

        #' @description
        #' every Reactome pathway has multiple reactions, each with a name
        #' @return list of names
       getReactionNames = function(){
          xml_text(xml_find_all(private$doc, "//reaction/@name"))
          },

        #' @description
        #' many pathways
        #' @return the map
       getPathwayNameMap = function(){
          private$tbl.pathwayNameMap
          },

        #' @description
        #' molecular species from xml hierarchy as list list
        #' @returns a named list, indexed by species ids, with reactome data in each element
       getMolecularSpeciesMap = function(){
         private$tbl.molecularSpecies
         },

        #' @description
        #' read the specified reaction, convert to edge and node tables
        #' @returns a named list, edges and nodes, each a data.frame
       processReaction = function(i, excludeUbiquitousSpecies, includeComplexMembers){
         reaction <- xml_find_all(private$doc, sprintf("//reaction[%d]", i))
         parser <- ReactionParser$new(private$doc, reaction, self)
         parser$toEdgeAndNodeTables(excludeUbiquitousSpecies, includeComplexMembers)
         }


       ) # public

) # class




