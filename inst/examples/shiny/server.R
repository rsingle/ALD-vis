#devtools::install_github("vpaunic/ALD-Measure")
library(shiny)
library(asymLD)

# Function to check whether package is installed
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 

if (!is.installed("fields")){
  install.packages("fields")
}


script <- "
//The following are d3 scaling functions.
//We are using a linear scale (scale.linear())
//There are 4 values in both the domain and range representing
//the four different colors represented in the color gradient
var homzScale = d3.scale.linear()
.domain( [0,         .33,        .66,        1])
.range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])

var freqScale = d3.scale.linear()
.domain( [0,         .33,        .66,        1])
.range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])

//This is a quick function to take a two value array of max and min values
//and spit out a four element array of evenly spaced values between those two points.
//If this was in R it would be easier (and may be easier in JS, I just dont know how.)
var domainGen = function(ar){
//Javascript indexes starting at 0, so this is the first value, or min.
var start = ar[0]
var end = ar[1]
var range = end - start
//return the finished array.
return [start, start + (0.33*range), start + (0.66*range), end]
}

//Waits to execute the javascript code 200 miliseconds. This is needed because shiny loads the javascript before it loads the table
//and thus has nothing to color.
window.setInterval(function() {

//Reach in and grab the max and min values from the hidden output renderText()
//See lines ~247-251 for the location of these values.
var maxFreq = parseFloat(d3.select('#maxVal_f').text())
var minFreq = parseFloat(d3.select('#minVal_f').text())
freqColorRange = [0, maxFreq]
//modify the domain of the frequency color scale from earlier with these newly
//obtained max and min values
freqScale.domain(domainGen(freqColorRange))

//Now we repeat this with the homz values
var maxHomz = parseFloat(d3.select('#maxVal_h').text())
var minHomz = parseFloat(d3.select('#minVal_h').text())
homzColorRange = [0, 1]
homzScale.domain(domainGen(homzColorRange))

//Now we use d3 to select all ('selectAll') of the frequency values by
//specifically targeting the 4th column of the html table.
//Then it sets the css style (.style) background-color by first grabbing the
//particular cell's value (d3.select(this).text()) and then putting that value
//into the color scaling function we set up in the above code.
//Set the colors for the allele freq
d3.selectAll('#asf_table tbody tr td:nth-child(4)')
.style('background-color', function() {
//grab the cells value
var cellValue = d3.select(this).text();
//return for the background-color value the result of the value run through color scale.
return (freqScale(cellValue))
})

//we then repeat this for the homz values
d3.selectAll('#asf_table tbody tr td:nth-child(5)')
.style('background-color', function() {
var cellValue = d3.select(this).text();
return (homzScale(cellValue))
})
//this is the amount of miliseconds between runs of the function. 
}, 200);
"

shinyServer(function(input, output, session) {
  
  session$onFlushed(function() {
    session$sendCustomMessage(type='jsCode', list(value = script))
  }, once = FALSE)
  
  dataInput <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    data <- read.csv(inFile$datapath, header=TRUE, sep=input$sep)
    # remove haplotypes with freq == 0 
    rm.inds <- data[, dim(data)[2]] == 0
    data <- data[!rm.inds, ]  
  })

  dataInput2 <- reactive({
    data <- dataInput()
    data.type2 <- TRUE
    names.type2 <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
    check.names <- names.type2 %in% names(data)    
    if (sum(!check.names) > 0) data.type2 <- FALSE
    if (is.null(data))
      return(NULL)
    else {

    if (data.type2){
      data$locus <- paste(data$locus1, data$locus2, sep="-")
      min.hf.sum <- min( aggregate(data$haplo.freq, by=list(data$locus), FUN=sum)[,2] )
      if ( min.hf.sum < 1 - input$tol) {
        output$text_rawdata1 <- renderText({ paste("") })
        output$text_rawdata2 <- renderText({ paste("Sum of haplo.freqs (min over locus pairs) = ", min.hf.sum, "(adjust tolerance value)") })
        return(NULL)
      }
      else {
        output$text_rawdata1 <- renderText({ paste("Sum of haplo.freqs (min over locus pairs) = ", min.hf.sum) })
        output$text_rawdata2 <- renderText({ paste("") })
        return(data)        
      }       
    } else {
      sum.freqs <- sum(data[, dim(data)[2]])
      if ( sum.freqs < 1 - input$tol) {
        output$text_rawdata1 <- renderText({ paste("") })
        output$text_rawdata2 <- renderText({ paste("Sum of haplo.freqs = ", sum.freqs, "(adjust tolerance value)") })
        return(NULL)
      }
      else {
        output$text_rawdata1 <- renderText({ paste("Sum of haplo.freqs = ", sum.freqs) })
        output$text_rawdata2 <- renderText({ paste("") })
        return(data)        
      }
    } 
  }  
  })  
  
  # display data
  output$raw_data <- renderDataTable({
    dataInput2()
  })
   
  # prepare the data for plotting function
  plotData <- reactive({
    withProgress(message = "Computing ALD values...", value = 0, {  
    #Sys.sleep(2) #useful for debugging progress bar stuff. 
    data <- dataInput2()
    if (is.null(data))
      return(NULL)
    else {
      data.type2 <- TRUE
      names.type2 <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
      check.names <- names.type2 %in% names(data)    
      if (sum(!check.names) > 0) data.type2 <- FALSE
      if (data.type2){
        data$locus <- paste(data$locus1, data$locus2, sep="-")
        levs <- unique(data$locus)
        data$locus <- factor(data$locus, levels=levs) #preserve ordering in the dataset        
        out<- by(data, list(locus=data$locus), compute.ALD)
        rbind.out <- NULL
        for (i in 1:length(out)) rbind.out <- rbind(rbind.out,out[[i]])
        ald.allpairs <- rbind.out 
      } else {
       #pop.use <- names(data)[dim(data)[2]]
        nloci <- dim(data)[2]-1
        loci <- names(data)[1:nloci]
        ald.allpairs <- NULL
        for (i in 1:(nloci-1)){
          for (j in (i+1):nloci){
            bi.data <- get_bilocus_data(data, i, j)
            ald.allpairs <- rbind(ald.allpairs, compute.ALD(bi.data, tolerance=input$tol))      
          }
        }        
      }
      ald.allpairs
    }
    })
  })
  
  # display statistic definitions above plotData
  output$text_ALDdata1 <- renderText({ paste(pre(
    "F.1     Homozygosity (expected under HWP) for locus 1",
    "F.2     Homozygosity (expected under HWP) for locus 2",
    "F.1.2   Conditional homozygosity for locus1 given locus2",
    "F.2.1   Conditional homozygosity for locus2 given locus1",
    "ALD.1.2 Asymmetric LD for locus1 given locus2",
    "ALD.2.1 Asymmetric LD for locus2 given locus1"
    )) 
  })
  
  # display plotData
  output$plot_data <- renderDataTable({
    plotData()
  })

  output$choose_locus_pair <- renderUI({
    data <- plotData()
    data$locus <- paste(data$locus1, data$locus2, sep="-")
    loci <- data$locus
    selectInput("selected_pair", "Choose a locus pair:", as.list(loci)) 
  })  
  output$choose_locus <- renderUI({
    if (is.null(plotData()))
      return(NULL)
    else {
      loci <- tryCatch({
        unlist(strsplit(input$selected_pair,"-"))
      }, warning = function(w) {
        #nothing
      }, error = function(e) {
        #nothing
      }, finally = {
        #all good
      }                
                        )
      selectInput("selected_locus", "Choose the focal locus:", as.list(loci)) 
    }
  })
  
  # plot of asymetric LD      
  output$heatmap <- renderPlot({
    withProgress(message = "Rendering Plot...", value = 0, {
    data <- plotData()
    if (is.null(data))
      return(NULL)
    else {
      loci <- unique(c(as.character(data$locus1),as.character(data$locus2)))
      ld.matrix.plot.2vars(dat=data, 
                           ld.varnames=c("ALD.1.2","ALD.2.1"), map.order=loci,
                           ld.labnames=c("",""), bw=T, xlab.shift=1, ylab.shift=-0, 
                           values=input$values)
      title(sub=paste("Asymmetric LD\n row gene conditional on
      column gene"),font.sub=2,cex.sub=1.2)            
    }
    })
  })
  
########################################################################################################################
# asf table
#Calculate the table inside of a reactive function. 
  calcTable <- reactive({
    data <- dataInput2()
    if (is.null(data))
      return(NULL)
    else {
      data.type2 <- TRUE
      names.type2 <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
      check.names <- names.type2 %in% names(data)    
      if (sum(!check.names) > 0) data.type2 <- FALSE
      if (data.type2){
        #loci <- c( unique(as.character(data$locus1)), unique(as.character(data$locus2)) )
        #bi.data <- data[data$locus1==loci[1] & data$locus2==loci[2],]
        loci <- unlist(strsplit(input$selected_pair,"-"))
        bi.data <- data[data$locus1==loci[1] & data$locus2==loci[2],]
        bi.data$locus1 <- as.character(bi.data$locus1)
        bi.data$locus2 <- as.character(bi.data$locus2)
        table <- compute.AShomz(bi.data, sort.var=c("focal","allele.freq"), sort.asc=c(F,F), tolerance=input$tol)
      } else {
        loci <- unlist(strsplit(input$selected_pair,"-"))
        loci.no <- (1:length(names(data)))[names(data) %in% loci]
        bi.data <- get_bilocus_data(data, loci.no[1], loci.no[2])
        bi.data$locus1 <- as.character(bi.data$locus1)
        bi.data$locus2 <- as.character(bi.data$locus2)
        table <- compute.AShomz(bi.data, sort.var=c("focal","allele.freq"), sort.asc=c(F,F), tolerance=input$tol)
      }
      table[table$focal==input$selected_locus,]        
    }
  })
  
  maxFreq <<- 0
  output$asf_table <- renderDataTable({
    #Run the table function inside of the renderDataTable
    #Wrapped it in a try catch because it was spilling out an error about undefined columns before finishing. 
    tryCatch({
      calcTable() 
    }, warning = function(w) {
      #nothing
    }, error = function(e) {
      #nothing
    }, finally = {
      #all good
    })
  }, options = list(orderClasses = TRUE))
  #Write the max and min values to the html, these will later be made hidden with css, They are for the javascript to use for
  #calculating the max and min for color scales. 
  output$maxVal_f <- renderText({ paste(max(calcTable()$allele.freq))})
  output$minVal_f <- renderText({ paste(0) }) #renderText({ paste(min(calcTable()$allele.freq))})
  output$maxVal_h <- renderText({ paste(1) }) #renderText({ paste(max(calcTable()$as.homz))})
  output$minVal_h <- renderText({ paste(0) }) #renderText({ paste(min(calcTable()$as.homz))})

########################################################################################################################

  output$asf_display <- renderUI({
    list(
      tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode", function(message) { eval(message.value); });')))
      , dataTableOutput("asf_table")
    )
  })
  
########################################################################################################################
  
  # -------------------------- HELPER FUNCTIONS -------------------------------
  
  # ---------------------------------------------------------------------------
  # Get the frequency data for two loci from the frequency file and format it 
  # for the compute.ALD() function.
  # ---------------------------------------------------------------------------
  get_bilocus_data <- function(data, i, j){
    nloci <- dim(data)[2]-1
    loci <- names(data)[1:nloci]
    data.2 <- data[,c(i, j, nloci+1)]
    data.2$pair <- paste(data.2[, 1], data.2[, 2], sep='-')
    aggregate.freqs <- by(data.2, data.2$pair, function(x) sum(x[, 3]))
    allele_combos <- names(aggregate.freqs)
    haplo.freq <- as.numeric(aggregate.freqs)
    alleles <- t(data.frame(strsplit(allele_combos, split='-')))
    colnames(alleles) <- c('allele1', 'allele2')
    bi.data <- cbind(data.frame(haplo.freq, locus1=rep(loci[i], length(haplo.freq)), 
      locus2=rep(loci[j], length(haplo.freq))), alleles)
    row.names(bi.data) <- NULL
    return(bi.data)
  }
  
  # ---------------------------------------------------------------------------
  # Function for plotting the ALD heatmap
  # ---------------------------------------------------------------------------
  
  ld.matrix.plot.2vars <- function(dat, ld.varnames, map.order, 
    ld.labnames=ld.varnames, bw=FALSE, xlab.shift=1, ylab.shift=-5, values=FALSE, cex.ld=.9) {
    dat2 <- dat
    
    # swap order of locus1 & locus2 in cases where they are not listed in 
    # map.order with locus1 coming before locus2
    dat2 <- dat2[dat2$locus1 %in% map.order & dat2$locus2 %in% map.order,]
    dat2$locus1 <- factor(dat2$locus1, levels=map.order)
    dat2$locus2 <- factor(dat2$locus2, levels=map.order)
    dat2$swap <- rep(0,dim(dat2)[1]) 
    dat2$swap[as.numeric(dat2$locus1) > as.numeric(dat2$locus2)] <- 1
    dat2$swap.locus1 <- dat2$locus1
    dat2$locus1[dat2$swap==1] <- dat2$locus2[dat2$swap==1]
    dat2$locus2[dat2$swap==1] <- dat2$swap.locus1[dat2$swap==1]  
    
    # make local copies of the vars on the dataframe to use from now on
    pval.color1 <- dat2[,names(dat2)==ld.varnames[1]]
    pval.color2 <- dat2[,names(dat2)==ld.varnames[2]]
    locus1 <- dat2$locus1
    locus2 <- dat2$locus2
    
    n.loci <- length(map.order)
    loci <- 1:length(map.order)
    names(loci) <- map.order
    
    x <- 1:n.loci
    y <- 1:n.loci
    names(x) <- map.order
    names(y) <- map.order
    
    # set up a matrix w/ rownames & colnames = loci.ordered
    z <- outer(x, y, FUN="+") 
    
    for (i in 1:(n.loci-1)) 
    {
      for (j in i:(n.loci))
      {
        if (i == j) 
        {
          z[i,j] <- -1 # a value outside of the zlim range for the diagonal
        } else 
        {
          # ld.varnames[1] on the upper triangle
          if (sum(locus1==names(x)[i] & locus2==names(y)[j]) > 0)  
          {
            #1st entry from pval.color1
            z[i,j] <- pval.color1[locus1==names(x)[i] & locus2==names(y)[j]] 
          } else
          {
            z[i,j] <- -1
          }
          # ld.varnames[2] on the lower triangle
          if (sum(locus1==names(x)[i] & locus2==names(y)[j]) > 0)  
          {
            #2nd entry from pval.color2
            z[j,i] <- pval.color2[locus1==names(x)[i] & locus2==names(y)[j]] 
          } else
          {
            z[j,i] <- -1
          }
        }  
      }
    }
    z[n.loci,n.loci] <- -1 # take care of bottom right element
    
    zz <- NULL
    for (i in (n.loci):1) 
    {
      zz <- cbind(zz, z[i,]) # rotate matrix z by +90 degrees
    }
    
    # NB: can change zz, but not z, since z is the return value
    zz[zz==-1] <- 9 # code untyped locus pairs as 9 for coloring
    
    # plot the LD values
    z.lim <- c(0,1) #zlim: plot only considers values between these 2 values
    stat.breaks <- c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
    if (bw)  stat.colors <- rev(gray((0:(length(stat.breaks)-2))/(length(stat.breaks)-2)))
    if (!bw) stat.colors <- c("gray",rev(rainbow((length(stat.breaks)-1)))[-(length(stat.breaks)-10)])
    # zlim: only use colors in this range
    fields::image.plot(x,y,zz,xaxt="n",yaxt="n",axes=F,xlab="",ylab="",zlim=z.lim,
      horizontal=F,col=stat.colors,breaks=stat.breaks) 
    
    # Add a light grid with dashed lines
    # NB: abline() must be called before mtext() and axis()
    # abline(h=0.5:(max(y)+.5),v=0.5:(max(x)+.5),lty=2,col=gray(.3))
    abline(h=0.5:(max(y)+.5),v=0.5:(max(x)+.5),lty=1,col="black")
    
    
    # add 1st var name as a text string rotated 90 degrees on the right margin
    mtext(side=4,ld.labnames[1],srt=90,font=2,line=ylab.shift) 
    # add 2nd var name as a text string on the bottom margin
    mtext(side=1,ld.labnames[2], font=2,line=xlab.shift) 
    
    # Add axes with locus names  
    axis(3, at=x, labels=names(x), tick=F,las=2)
    axis(2, at=y, labels=rev(names(y)), tick=F,las=2) # rev() since rotated matrix 90 degrees 
    
    
    if (values)
    {
      for (x in 1:n.loci)
      {
        for (y in n.loci:1)
        {
          if (x != (n.loci-y+1))
          {
            yyz <- as.character(round(zz[x,y],2))
            if (!is.na(zz[x,y]))
            {
              if (nchar(as.character(yyz))==3) yyz <- paste(as.character(yyz),"0",sep="")
              yyz <- substr(yyz,2,4)
              if (zz[x,y]>.5 & bw==T)
                # font=2 for bold, NB: cex was 1.2
                text(x,y,yyz,cex=cex.ld,col="white",font=1) 
              else 
                # font=2 for bold, NB: cex was 1.2
                text(x,y,yyz,cex=cex.ld,col="black",font=1) 
            }
          }
        }
      }
    }
    return(z)
  }
  
  
  })

