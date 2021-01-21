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
              elements.raw <- xml_text(xml_find_all(species, ".//bqbiol:is//rdf:li/@rdf:resource"))
              elements <- c()
              chebi.id <- NA
              uniprot.id <- NA
              if(length(elements.raw) > 0){
                  chebi.element <- grep("CHEBI", elements.raw)
                  uniprot.element <- grep("uniprot", elements.raw)
                  if(length(chebi.element) == 1){
                      chebi.id <- sub("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", "",
                                      elements.raw[chebi.element], fixed=TRUE)
                      moleculeType <- "smallMolecule"
                      }
                  if(length(uniprot.element) == 1){
                      uniprot.id <- sub("http://purl.uniprot.org/uniprot/", "",
                                        elements.raw[uniprot.element], fixed=TRUE)
                      moleculeType <- "protein"
                      }
                 } # if some elements
              members <- xml_text(xml_find_all(species, ".//bqbiol:hasPart//rdf:li/@rdf:resource"))
              is.elements <- xml_text(xml_find_all(species, ".//bqbiol:is//rdf:li/@rdf:resource"))
              notes <- xml_text(xml_find_all(species, ".//notes"))
                # notes are sometimes missing.  assume species is a protein.  might be a terrible idea
              if(length(notes) == 0){
                 printf("=============================================")
                 printf("no notes section for %s in %s", name, private$sbml.filename)
                 printf("=============================================")
                 alternativeEntities <- FALSE
                 complex <- FALSE
                 protein <- TRUE
                 smallCompound <- FALSE
                 }
              if(length(notes) > 0){
                 alternativeEntities <- grepl("alternative entities", notes)
                 complex <- grepl("Reactome Complex", notes)
                 protein <- grepl("This is a protein", notes)
                 smallCompound <- grepl("This is a small compound", notes)
                }

              printf("--- name: %s", name)
              #if(name == "SLC36A4") browser()
              #print(protein)
              #print(members)
              #print(is.elements)
              if(protein & length(members) == 0 & length(is.elements) > 0){
                 moleculeType <- "protein"
                 uniprot.id <- sub("http://purl.uniprot.org/uniprot/", "uniprotkb:", is.elements[1])
                 members <- uniprot.id
                 }
              #if(length(members) == 0)
              #  members <- c()
                # a pseudo complex, where the plural members are actually alternate proteins?
                # todo: need a more robust solution here. now just use the first alternate
              #printf("length(members): %d   ; found? %s", length(members), !grepl(";", name))
              if(length(members) > 0 & alternativeEntities){
                  moleculeType <- "protein"
                  uniprot.id <- sub(".*http://purl.uniprot.org/uniprot/", "uniprotkb:", members)[1]
                  uniprot.id <- sub('"', '', uniprot.id)
                  members <- uniprot.id
                  }
              if(length(members) > 0 & complex){
                moleculeType <- "complex"
                members <- sub("http://purl.uniprot.org/uniprot/", "uniprotkb:", members, fixed=TRUE)
                members <- sub("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:", "ChEBI:", members, fixed=TRUE)
                members <- sub("http://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=",
                               "ligandId:", members, fixed=TRUE)
                }
              new.entry <- list(name=name, moleculeType=moleculeType, compartment=compartment,
                                members=members, uniprot=uniprot.id, chebi=chebi.id)
              tmp[[id]] <- new.entry
              } # for species
           extract <- function(i){
               x <- tmp[[i]]
               complex.members <- x$members
               data.table(id=names(tmp)[i], name=x$name, type=x$moleculeType,
                          compartment=x$compartment,
                          uniprot.id=x$uniprot,
                          chebi.id=x$chebi,
                          complex.members=list(complex.members))
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
        #' molecular species from xml hierarchy as list
        #' @returns a named list, indexed by species ids, with reactome data in each element
       getMolecularSpeciesMap = function(){
         private$tbl.molecularSpecies
         },

        #' @description
        #' read the specified reaction, convert to edge and node tables
        #' @param i integer, which reaction from the current pathway
        #' @param excludeUbiquitousSpecies logical, e.g. water, atp, adp
        #' @param includeComplexMembers logical, expand graph to include molecules participating in complexes
        #' @returns a named list, edges and nodes, each a data.frame
       processReaction = function(i, excludeUbiquitousSpecies, includeComplexMembers){
         reaction <- xml_find_all(private$doc, sprintf("//reaction[%d]", i))
         parser <- ReactionParser$new(private$doc, reaction, self)
         parser$toEdgeAndNodeTables(excludeUbiquitousSpecies, includeComplexMembers)
         }


       ) # public

) # class




