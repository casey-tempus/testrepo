

library(RMySQL)
mydb = dbConnect(MySQL(), user='dbreader', password='dbreader', dbname='outside_db', host='172.22.10.42')
dbListTables(mydb)
dbListFields(mydb, 'cosmic')
# SELECT TOP 10 distinct * 
#  FROM people 
# WHERE names='SMITH'
# ORDER BY names asc

rs = dbSendQuery(mydb, "SELECT * FROM cosmic LIMIT 10")
data = fetch(rs, n=-1)
Gene
Mutation_CDS
Mutation_AA

pubs  <- pubs.clinical("GNAS", list("NM_000516", "c.602G>A", "p.Arg201His", "20:57484421G>A"))
  
library(data.table)
# function for searching pubmed directly from keywords  
pubs.clinical        <- function ( gene, annotations ) {
  
    pubs.var            <- do.call("rbind", 
                                   lapply( 1:length( annotations ), 
                                           function(i) get.pubmed.ids( c( gene,  annotations[ i ] ))))
    
    pubs.var.freq       <- sort( table( pubs.var$pmid ), decreasing = T )
    pubs.var.sorted     <- pubs.var[ order( match( pubs.var$pmid, names( pubs.var.freq )), decreasing = F ), ]
    pubs.var.summary    <- setDT( pubs.var.sorted )[ , lapply( .SD, function(x) paste( rev( unique( unlist( strsplit( gsub( "^\\s+|\\s+$", "", toString( na.omit(x) )), ",")))), collapse =" + ")), by = list( pmid, date, cited, title, journal )]
    results <- pubs.var.summary
    
    return(results)
  
}



library(RISmed)
# function for searching pubmed directly from keywords
get.pubmed.ids         <- function(search.elements){
  
  require(RISmed)
  
  startDate             <- 2000
  endDate               <- 2017
  ptm                   <- proc.time()
  search.phrase         <- paste( search.elements, sep = "", collapse = " AND " )
  search.query          <- EUtilsSummary( search.phrase, mindate=startDate, maxdate=endDate )
  records               <- EUtilsGet( search.query )
  cites                 <- Cited( records )
  dates                 <- YearPubmed( records )
  journals              <- ISOAbbreviation( records )
  titles                <- ArticleTitle( records )
  pmids                 <- names( cites )
  
  #author afilliation
  #print(Affiliation(records))
  #print(MedlineTA(records))
  #print(ISOAbbreviation(records))
  
  search.comb           <- paste( search.elements, collapse=' , ' )
  search                <- rep( search.comb, length(pmids))
  
  pmids.sort            <- pmids[    order( cites, decreasing = TRUE )]
  dates.sort            <- dates[    order( cites, decreasing = TRUE )]
  cites.sort            <- cites[    order( cites, decreasing = TRUE )]
  titles.sort           <- titles[   order( cites, decreasing = TRUE )]
  journals.sort         <- journals[ order( cites, decreasing = TRUE )]
  
  results               <- data.frame( 'search'  = search, 
                                       'pmid'    = pmids.sort, 
                                       'date'    = dates.sort , 
                                       'cited'   = cites.sort,
                                       'journal' = journals.sort,
                                       'title'   = titles.sort )
  
  return( results )
  
}