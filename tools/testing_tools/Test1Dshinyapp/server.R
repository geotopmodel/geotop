#
#
#
#
rm(list=ls())

library(geotopbricks)
library(dygraphs)


inpts.file <- "geotop.inpts"
##wpath_1D <- '/home/ecor/local/src/geotop_dev/tests/1D' 
wpath_1D <- '../../../tests/1D' 
wpaths <- list.files(wpath_1D,recursive=FALSE,full.name=TRUE)
names(wpaths) <- list.files(wpath_1D,recursive=FALSE,full.name=FALSE)

wpaths <- wpaths[names(wpaths)=="Matsch_B2_Ref_007"]

## Check Simulation Occurence 

inpts.files <- paste(wpaths,inpts.file,sep="/")
igeotopsim <- which(file.exists(inpts.files))

wpaths <- wpaths[igeotopsim]
inpts.files <- inpts.files[igeotopsim]


#wpaths <- wpaths[names(keywords.out)!="Calabria"]
#wpaths <- wpaths[names(keywords.out)!="ColdelaPorte"]
#wpaths <- wpaths[names(keywords.out)!="Jungfraujoch"]



#### Search the output
output_values <- c("output-tabs") ###,"output-maps")

keywords <- lapply(X=wpaths,FUN=declared.geotop.inpts.keywords,inpts.file=inpts.file)

timezone <- sapply(X=keywords,FUN=function(x){
			
	        o1 <- str_trim(x$Value[x$Keyword=="MeteoStationStandardTime"])[1]
			o2 <- str_trim(x$Value[x$Keyword=="StandardTimeSimulation"])[1]
			
			print(o1)
			print(o2)
			o <- c(o1,o2)
			o <- o[!is.na(o)][1]
			o <- -as.integer(o)
			if (o>=0) {
				o <- sprintf("Etc/GMT+%d",o)
			} else {
				
				o <- sprintf("Etc/GMT%d",o)
			}
			return(o)
			
		})


keywords.out <- lapply(X=keywords,FUN=function(x,keyword_values){
			cnt <- NULL
			for (kit in keyword_values) {
				
				cnt <- union(cnt,which(str_detect(x$Value,kit)))
				
			}
			
			x[cnt,]
			},keyword_values=output_values)


for (it in names(keywords.out)) {
	
	
	keywords.out[[it]]$sim_name <- it
	keywords.out[[it]]$sim_wpath <- wpaths[it]
	keywords.out[[it]]$formatter <- "%04d"
	keywords.out[[it]]$formatter[keywords.out[[it]]$Keyword %in% c("BasinOutputFile","DischargeFile","SnowCoveredAreaFile")] <- ""
	keywords.out[[it]]$tz <- timezone[[it]]
	keywords.out[[it]]$level <- 1 ##paste(c(1,32,33),collapse=",")
	if (it=="Jungfraujoch") keywords.out[[it]]$level <- "32,33"
}

#keywords.out <- keywords.out[names(keywords.out)!="Calabria"]
#keywords.out <- keywords.out[names(keywords.out)!="ColdelaPorte"]
#keywords.out <- keywords.out[names(keywords.out)!="Jungfraujoch"]
#stop("HERE")
#####
# suffixes <- c("-SE27XX","-METEOIO-ON","-METEOIO-OFF","")
suffixes <- c("-v_2.0","-METEOIO-ON","-v_3.0","")
names(suffixes) <- suffixes
names(suffixes)[suffixes==""] <- "latest"
values.out <- list() 

#print("ccc")
#print(suffixes)
#print(keywords.out)

for (it_s in names(suffixes)) {
#	print(it_s)
	values.out[[it_s]] <- lapply(X=keywords.out,FUN=function(x,add_suffix_dir,inpts.file) {

		#	print(x)
			vars <- x$Keyword
			wpath <- x$sim_wpath[1]
			formatter <- x$formatter
			tz <- x$tz
			
			level <- x$level
	#		print(level)
			if (is.character(level)) {
				level <- str_split(level,",")
				print(level)
				
			} else {
				
				level <- as.list(level)
			}
			
	
			o <- list()
			
			
			for (i in 1:length(vars)) {
				
				lic <- 1:length(level[[i]])
				if (formatter[i]=="") lic <- 1
				
				date_field <- "Date12.DDMMYYYYhhmm."
				if (vars[i]=="DischargeFile") date_field <- "DATE.day.month.year.hour.min."
				for (li in lic) {
				lll <- as.integer(level[[i]][li])
				itn <- sprintf("%s%04d",vars[i],lll)
				
				if (formatter[i]=="") itn <- vars[i]
				o[[itn]] <- list(keyword=vars[i],inpts.file=inpts.file,
				wpath=wpath,data.frame=TRUE,date_field=date_field,tz=tz[1],
				formatter=formatter[i],add_suffix_dir=add_suffix_dir,level=lll)
				v <- o[[itn]]
				v[["header.only"]] <- TRUE
				attr(o[[itn]],"header") <- try(do.call(what=get.geotop.inpts.keyword.value,args=v),silent=TRUE)

				}
			}
			
			oc <- sapply(X=o,FUN=function(x) {class(attr(x,"header"))})	
			
			o <- o[which(oc!="try-error")] 
			return(o)
			
			
		},add_suffix_dir=suffixes[it_s],inpts.file=inpts.file)


}

sim_names <- names(values.out[[1]])

# ...


wpathm <- sapply(X=keywords.out,FUN=function(x){x$sim_wpath[1]})

#### METEO DATA
meteodata <- list()
for (itme in sim_names) {
	level <- 1
	wpme <- keywords.out[[itme]]$sim_wpath[1]
	tz <- keywords.out[[itme]]$tz[1]
	date_field <- try(get.geotop.inpts.keyword.value("HeaderDateDDMMYYYYhhmmMeteo",wpath=wpme,inpts.file=inpts.file,exceptions="none"),silent=TRUE)
	if (class(date_field)=="try-error") date_field <- "Date"
	
	start_date <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpme,tz=tz) 
	end_date <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpme,tz=tz) 
	
	
	meteo.args  <- list(keyword="MeteoFile",wpath=wpme,data.frame=TRUE,level=level,date_field=date_field,tz=tz,inpts.file=inpts.file,start_date=start_date,end_date=end_date,MAXNROW=1)
	meteodata[[itme]]  <- meteo.args   
	meteo.args[["header.only"]] <- TRUE
	
	attr(meteodata[[itme]],"header") <- try(do.call(what=get.geotop.inpts.keyword.value,args=meteo.args),silent=TRUE)
	
	if (class(attr(meteodata[[itme]],"header"))=="try-error")    {
		
		msg <- sprintf("Issue in meteo dat%s",paste(unlist(meteodata[[itme]]),collapse=" "))
		stop(msg)
		
	}                                                       
		
	
	
	
	
}



## SERVER


shinyServer(function(input, output) {
#	
#	output$scatterPlot <- renderPlot({
#				x <- rnorm(input$n)
#				y <- rnorm(input$n)
#				plot(x, y)
#			})
	##output$sim_vars <- reactive(c("LOADING",sim_vars))
##	output$sim_vars_ui  <-renderUI({selectInput('var', 'Variable',c("LOADING",sim_vars))})

  
	output$sim_names_ui <-renderUI({selectInput('sim', 'GEOtop Simulation Test',c("LOADING",sim_names))})
	output$sim_kws_ui <-renderUI({selectInput('kws', 'Variable / Keyword',c("LOADING",names(values.out[[1]][[input$sim]])))})
	output$sim_layer_ui <-renderUI({selectInput('layer', 'Layer / Variable',c("LOADING",attr(values.out[[1]][[input$sim]][[input$kws]],"header")))})
	output$sim_meteo_ui <-renderUI({selectInput('meteo', 'Meteorological Forcing Variable',c("LOADING",attr(meteodata[[input$sim]],"header")))})
	###uiOutput("dupes")
	
	
	output$dygraph <- 
			renderDygraph({
						
					suffixes_ <- names(suffixes)
					if (input$latest==FALSE) suffixes_ <- suffixes_[suffixes_!="latest"]
					if (input$SE27XX==FALSE) suffixes_ <- suffixes_[suffixes_!="-v_2.0"]
					if (input$METEOIO_ON==FALSE) suffixes_ <- suffixes_[suffixes_!="-METEOIO-ON"]
					if (input$METEOIO_OFF==FALSE) suffixes_ <- suffixes_[suffixes_!="-v_3.0"]
					
					
					icsf <- which(names(values.out) %in% suffixes_)	
						
					sim <- input$sim 
					kws <- input$kws
					layer <- input$layer
					data_l <- lapply(X=values.out[icsf],FUN=function(x,sim,kws,layer){
								o <- try(do.call(args=x[[sim]][[kws]],what=get.geotop.inpts.keyword.value),silent=TRUE)
								if (class(o)=="try-error")  { 
									o <- NULL
								} else { 
								#print(class(o))
								#print(layer)
									o <- o[,layer]
								}
								#str(o)
								return(o)
								},
								sim=sim,kws=kws,layer=layer)
					print(data_l)
					
					data_time <- index(data_l[[1]])
					names_n <- names(data_l)
					data_l <- lapply(X=data_l,FUN=as.vector)
					str(data_l)
					data <- as.zoo(as.data.frame(do.call(what="cbind",args=data_l)))
					index(data) <- data_time 
						
					###	residuals <- data-data[,"-SE27XX"]
						
					str(data)
					str(residuals)
					colors <- RColorBrewer::brewer.pal(4, "Set1")[icsf]
						
					main <- paste(sim,kws,layer,sep="::")	
					
					dygraph(data, ylab="[unit]",main=main) %>% dyRangeSelector() %>%
								dyRoller() %>%
								dyOptions(colors = colors,digitsAfterDecimal=input$digitsAfterDecimal) %>%
								dyLegend(labelsSeparateLines = TRUE)
						
					})
					
					
			output$dygraph_res <-	
					
					renderDygraph({
								
								suffixes_ <- names(suffixes)
								if (input$latest==FALSE) suffixes_ <- suffixes_[suffixes_!="latest"]
								if (input$SE27XX==FALSE) suffixes_ <- suffixes_[suffixes_!="-v_2.0"]
								if (input$METEOIO_ON==FALSE) suffixes_ <- suffixes_[suffixes_!="-METEOIO-ON"]
								if (input$METEOIO_OFF==FALSE) suffixes_ <- suffixes_[suffixes_!="-v_3.0"]
								
								
								icsf <- which(names(values.out) %in% suffixes_)	
								
								sim <- input$sim 
								kws <- input$kws
								layer <- input$layer
								data_l <- lapply(X=values.out[icsf],FUN=function(x,sim,kws,layer){
											o <- try(do.call(args=x[[sim]][[kws]],what=get.geotop.inpts.keyword.value),silent=TRUE)
											if (class(o)=="try-error")  { 
												o <- NULL
											} else { 
												#print(class(o))
												#print(layer)
												o <- o[,layer]
											}
											#str(o)
											return(o)
										},
										sim=sim,kws=kws,layer=layer)
								print(data_l)
								
								data_time <- index(data_l[[1]])
								names_n <- names(data_l)
								data_l <- lapply(X=data_l,FUN=as.vector)
								str(data_l)
								data <- as.zoo(as.data.frame(do.call(what="cbind",args=data_l)))
								index(data) <- data_time 
								
								###	residuals <- data-data[,"-SE27XX"]
								
								#str(data)
								#str(residuals)
								colors <- RColorBrewer::brewer.pal(4, "Set1")[icsf]
								
								residuals <- data-data[,"-v_2.0"]
								
								main <- paste(sim,kws,layer,sep="::")	
								dygraph(residuals, ylab="[unit]",main=main) %>% dyRangeSelector() %>%
										dyRoller() %>%
										dyOptions(colors = colors,digitsAfterDecimal=input$digitsAfterDecimal) %>%
										dyLegend(labelsSeparateLines = TRUE)
								
								
							})
					
							
		output$dygraph_meteo <-	
									
				renderDygraph({
							    meteocol <- input$meteo
							    sim <- input$sim 
								data <- try(do.call(args=meteodata[[sim]],what=get.geotop.inpts.keyword.value),silent=TRUE)
								if (class(data)=="try-error")  { 
									data <- NULL
								} else { 
									#print(class(o))
									#print(layer)
									data <- data[,meteocol]
								}	
							
								data[data<=-9990.0] <- NA
								main <- paste(sim,"meteo0001",meteocol,sep="::")
								dygraph(data, ylab="[unit]",main=main) %>% dyRangeSelector() %>%
										dyRoller() %>%
										dyOptions(digitsAfterDecimal=input$digitsAfterDecimal) %>%
										dyLegend(labelsSeparateLines = TRUE)
												
												
												
				})
							
					#renderDygraph({
								
						
								
							
								
					#		})

			  output$info <- renderText({
								"bla-bla"
							})
			  
	
})

#
#
#shinyApp(
#		ui = fluidPage(
#				fluidRow(
#					selectInput('var', 'Variable',sim_vars)
#					#	column(4, selectInput('ycol', 'Y Variable', names(dataset),
#					#					selected=names(dataset)[[2]])),
#					#	column(4, numericInput('clusters', 'Cluster count', 3,
#					##					min = 1, max = 9))
#				),
#				fluidRow(
#						plotOutput('kmeans', height = "400px")  
#				)
#		),
#		
#		server = function(input, output, session) {
#			
#			# Combine the selected variables into a new data frame
##			selectedData <- reactive({
##						dataset[, c(input$xcol, input$ycol)]
##					})
##			
##			clusters <- reactive({
##						kmeans(selectedData(), input$clusters)
##					})
##			
##			output$kmeans <- renderPlot(height = 400, {
##						par(mar = c(5.1, 4.1, 0, 1))
##						plot(selectedData(),
##								col = clusters()$cluster,
##								pch = 20, cex = 3)
##						points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
##					})
#		},
#		
#		options = list(height = 500)
#)
#
#stop("CIAO")
#
#

#
#
#
#
#
#
#kmeans_cluster <- function(dataset) { 
#	
#	require(shiny)  
#	
#	shinyApp(
#			ui = fluidPage(
#					fluidRow(style = "padding-bottom: 20px;",
#							column(4, selectInput('xcol', 'X Variable', names(dataset))),
#							column(4, selectInput('ycol', 'Y Variable', names(dataset),
#											selected=names(dataset)[[2]])),
#							column(4, numericInput('clusters', 'Cluster count', 3,
#											min = 1, max = 9))
#					),
#					fluidRow(
#							#plotOutput('kmeans', height = "400px")  
#					)
#			),
#			
#			server = function(input, output, session) {
#				
##				# Combine the selected variables into a new data frame
##				selectedData <- reactive({
##							dataset[, c(input$xcol, input$ycol)]
##						})
##				
##				clusters <- reactive({
##							kmeans(selectedData(), input$clusters)
##						})
##				
##				output$kmeans <- renderPlot(height = 400, {
##							par(mar = c(5.1, 4.1, 0, 1))
##							plot(selectedData(),
##									col = clusters()$cluster,
##									pch = 20, cex = 3)
##							points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
##						})
#			},
#			
#			options = list(height = 500)
#	)
#}
#
#
#
#shinyApp(
#		ui = fluidPage(
#				shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
#				sidebarPanel(
#						textInput("txt", "Text input:", "text here"),
#						selectInput(inputId = "variable", label = "discover variable", choices = unique(c("mela","pera"))),
#						sliderInput("slider", "Slider input:", 1, 100, 30),
#						actionButton("action", "Button"),
#						actionButton("action2", "Button2", class = "btn-primary")
#				),
#				mainPanel(
#						tabsetPanel(
#								tabPanel("Tab 1"),
#								tabPanel("Tab 2")
#						)
#				)
#		),
#		server = function(input, output) {
#			
#			#output$variable <- unique(c(input$variable,input$txt))
#			####ss <<- output$vaqriable
#			####print(output)
#		}
#)
#
#
#
######
##keywords.out <- melt(keywords.out)
##
#
#
#
#
#
#
#
#
